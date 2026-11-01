"""Strict HEST benchmark reader. Never assumes patch and expression row order."""
from dataclasses import dataclass
from pathlib import Path
import h5py
import numpy as np

@dataclass
class Slide:
    sample_id: str
    barcodes: np.ndarray
    genes: np.ndarray
    counts: np.ndarray
    xy_px: np.ndarray
    xy_um: np.ndarray
    thumbnail: np.ndarray
    thumbnail_scale: float
    patch_path: Path
    excluded_expression_locations: int
    coordinate_error_px: float
    excluded_control_features: list[str]


def strings(dataset):
    return np.asarray(dataset[:]).reshape(-1).astype(str)


def read_slide(root, sample_id, microns_per_pixel):
    if not np.isfinite(microns_per_pixel) or microns_per_pixel <= 0:
        raise ValueError('A positive calibrated microns-per-pixel value is required')
    patch_path = Path(root) / 'IDC/patches' / f'{sample_id}.h5'
    with h5py.File(Path(root) / 'IDC/adata' / f'{sample_id}.h5ad') as a, h5py.File(patch_path) as p:
        expr_ids, patch_ids = strings(a['obs/_index']), strings(p['barcode'])
        if len(set(expr_ids)) != len(expr_ids) or len(set(patch_ids)) != len(patch_ids):
            raise ValueError('Duplicate barcodes')
        lookup = {b: i for i, b in enumerate(expr_ids)}
        if not set(patch_ids) <= set(expr_ids):
            raise ValueError('Patch barcode missing from expression')
        indices = np.array([lookup[b] for b in patch_ids])
        xy = a['obsm/spatial'][:][indices].astype(float)
        coords = p['coords'][:]
        if coords.shape != xy.shape or not np.isfinite(xy).all():
            raise ValueError('Invalid coordinate shape or values')
        # HEST IDC patch coords are (row, column); AnnData spatial is (x, y).
        error = float(np.max(np.abs(coords[:, ::-1] - xy)))
        if error > 1.01:
            raise ValueError(f'Patch/expression coordinate mismatch: {error} pixels')
        counts = a['X'][:][indices]
        if not np.isfinite(counts).all() or np.any(counts < 0) or not np.allclose(counts, np.rint(counts)):
            raise ValueError('Expected nonnegative integer raw counts in X')
        genes = strings(a['var/_index'])
        if len(set(genes)) != len(genes):
            raise ValueError('Duplicate gene identifiers')
        base = a['uns/spatial/ST']
        thumb = base['images/downscaled_fullres'][:]
        scale = float(base['scalefactors/tissue_downscaled_fullres_scalef'][()])
        if len(p['img']) != len(patch_ids) or counts.shape != (len(patch_ids), len(genes)):
            raise ValueError('Inconsistent row dimensions')
        keep = np.array([not g.startswith(('NegControl', 'BLANK_', 'Unassigned', 'Deprecated')) for g in genes])
        excluded = genes[~keep].tolist()
        return Slide(sample_id, patch_ids, genes[keep], counts[:, keep], xy, xy * microns_per_pixel,
                     thumb, scale, patch_path, len(expr_ids)-len(patch_ids), error, excluded)


def validate_patient_splits(splits, patients):
    seen_slides, seen_patients = set(), set()
    if set(splits) != {'train', 'validation', 'test'}:
        raise ValueError('Expected train, validation and test splits')
    for name, slides in splits.items():
        if not slides or len(set(slides)) != len(slides):
            raise ValueError(f'Empty or duplicate slide list: {name}')
        current = {patients[s] for s in slides}
        if seen_slides.intersection(slides) or seen_patients.intersection(current):
            raise ValueError('Patient or slide leakage across splits')
        seen_slides.update(slides)
        seen_patients.update(current)
    if seen_slides != set(patients):
        raise ValueError('Split does not cover cohort exactly')


def harmonize_genes(slides):
    """Intersect assay identifiers only; never use held-out expression values.

    Name-based harmonization is essential: the HEST IDC files have different
    control sets and biological panel sizes despite equal total feature counts.
    """
    if not slides:
        raise ValueError('No slides to harmonize')
    shared = sorted(set.intersection(*(set(s.genes) for s in slides.values())))
    if not shared:
        raise ValueError('No shared biological genes')
    records = {}
    for sample, slide in slides.items():
        lookup = {g:i for i,g in enumerate(slide.genes)}
        records[sample] = {'biological_features_before_harmonization':len(slide.genes),
                           'nonshared_biological_features':sorted(set(slide.genes)-set(shared))}
        slide.counts = slide.counts[:,[lookup[g] for g in shared]]
        slide.genes = np.array(shared)
    return records
