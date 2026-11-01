"""Freeze C0 artifacts. Exits nonzero if the scientific audit remains incomplete."""
import hashlib
import json
import platform
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import pandas as pd
from graphpathology.data import read_slide, validate_patient_splits, harmonize_genes
from graphpathology.masks import REGIMES, make_mask, stable_seed, anchor_distances
from graphpathology.preprocessing import CountTransform, select_panels


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def main():
    start = time.monotonic()
    out = ROOT/'reports/C0'
    if (out/'C0_SPEC_LOCK.json').exists():
        raise SystemExit('C0 is frozen. Use a new versioned checkpoint directory for changes.')
    config = json.loads((ROOT/'configs/c0.json').read_text())
    metadata = pd.read_csv(ROOT/'data/raw/HEST_v1_1_0.csv').set_index('id')
    ids = sum(config['splits'].values(), [])
    patients = {s:str(metadata.loc[s,'patient']) for s in ids}
    validate_patient_splits(config['splits'], patients)
    write_json(out/'C0_SPLITS.json', {'splits':config['splits'], 'patients':patients,
                                    'patient_source':'Pinned HEST_v1_1_0.csv; global HEST patient IDs',
                                    'limitations':'Only four patients. Two train, one validation, one test; insufficient for publication-level inference.'})
    slides = {s:read_slide(ROOT/'data/raw',s,float(metadata.loc[s,'pixel_size_um_estimated'])) for s in ids}
    harmonization = harmonize_genes(slides)
    genes = slides[ids[0]].genes
    for slide in slides.values():
        if not np.array_equal(slide.genes, genes):
            raise ValueError('Gene sets/order differ; explicit harmonization required')
    train = np.concatenate([slides[s].counts for s in config['splits']['train']])
    transform = CountTransform().fit(train)
    panels = select_panels(train,genes,config['panel_sizes'])
    write_json(out/'C0_GENE_PANELS.json', {'selection':'Training-only variance of log1p counts; no formal HVG/SVG claim', 'panels':panels})
    write_json(out/'C0_NORMALIZATION.json', {'method':config['normalization'], 'genes':genes.tolist(),
               'mean':transform.mean.tolist(), 'scale':transform.scale.tolist(), 'fit_slides':config['splits']['train']})
    mask_dir = out/'C0_MASKS'; mask_dir.mkdir(exist_ok=True)
    mask_records, slide_records = [], []
    for sample, slide in slides.items():
        feat_path = ROOT/'data/cache'/f'{sample}_features.npz'
        with np.load(feat_path,allow_pickle=False) as f:
            if not np.array_equal(f['barcodes'],slide.barcodes):
                raise ValueError('Feature cache barcode order mismatch')
            for key in ['embeddings','morphology']:
                if len(f[key]) != len(slide.counts) or not np.isfinite(f[key]).all():
                    raise ValueError('Invalid feature cache')
            dims = {k:list(f[k].shape) for k in ['embeddings','morphology']}
        normalized = ROOT/'data/cache'/f'{sample}_targets.npz'
        np.savez_compressed(normalized, raw_counts=slide.counts,
                            normalized_targets=transform.transform(slide.counts),
                            barcodes=slide.barcodes, genes=genes, xy_um=slide.xy_um)
        for regime in REGIMES:
            for q in config['coverage_fractions']:
                for rep in range(config['mask_replicates']):
                    seed = stable_seed(sample,regime,q,rep,config['seed'])
                    mask, record = make_mask(slide.xy_um,regime,q,seed)
                    name = f'{sample}_{regime}_q{q:.2f}_r{rep}'
                    path = mask_dir/f'{name}.npz'
                    np.savez_compressed(path,observed=mask,barcodes=slide.barcodes,
                                        distance_to_anchor_um=anchor_distances(slide.xy_um,mask))
                    exposed, _ = transform.observed_inputs(slide.counts, mask)
                    corrupt = slide.counts.astype(float).copy(); corrupt[~mask] = np.nan
                    alternate, _ = transform.observed_inputs(corrupt, mask)
                    if not np.array_equal(exposed, alternate) or np.any(exposed[~mask]):
                        raise AssertionError('Hidden target leakage')
                    record.update(mask_id=name, sample_id=sample, replicate=rep,
                                  path=str(path.relative_to(ROOT)),sha256=sha(path),
                                  hidden_target_invariance=True)
                    mask_records.append(record)
        slide_records.append({'sample_id':sample,'patient':patients[sample],
            'split':next(k for k,v in config['splits'].items() if sample in v),
            'n_locations':len(slide.counts),'n_genes':len(genes),
            'excluded_control_features':slide.excluded_control_features,
            'gene_harmonization':harmonization[sample],
            'expression_locations_without_patches':slide.excluded_expression_locations,
            'zero_count_locations':int((slide.counts.sum(axis=1)==0).sum()),
            'max_coordinate_error_px':slide.coordinate_error_px,
            'microns_per_pixel':float(metadata.loc[sample,'pixel_size_um_estimated']),
            'physical_calibration':'HEST estimated pixel size; uncertainty not quantified',
            'feature_cache':str(feat_path.relative_to(ROOT)),'feature_sha256':sha(feat_path),
            'feature_shapes':dims,'feature_provenance':json.loads(feat_path.with_suffix('.json').read_text()),
            'target_cache':str(normalized.relative_to(ROOT)),'target_sha256':sha(normalized),
            'source_study':str(metadata.loc[sample,'study_link']).strip(),
            'source_page':str(metadata.loc[sample,'download_page_link1']).strip()})
    downloads = json.loads((out/'C0_DOWNLOADS.json').read_text())
    for record in downloads['files']:
        if sha(ROOT/record['path']) != record['sha256']:
            raise ValueError('Raw asset checksum mismatch')
    write_json(out/'C0_MASK_INDEX.json',mask_records)
    write_json(out/'C0_DATA_MANIFEST.json',{'dataset':config['dataset'], 'downloads':downloads,
                                          'slides':slide_records,'masks':len(mask_records)})
    pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in slide_records]).to_csv(out/'C0_DATASET_STATISTICS.csv',index=False)
    audit_path = out/'C0_ALIGNMENT_REVIEW.json'
    audit = json.loads(audit_path.read_text()) if audit_path.exists() else {'status':'pending'}
    passed = audit.get('status') == 'passed_thumbnail_and_patch_audit'
    source_commit = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    dirty = bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip())
    source_hashes = {str(p.relative_to(ROOT)):sha(p) for folder in ['src','scripts','configs','tests']
                     for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p)}
    spec = {'checkpoint':'C0', 'status':'PASS_DEVELOPMENT_ONLY' if passed else 'BLOCKED_ALIGNMENT_REVIEW',
            'timestamp':datetime.now(timezone.utc).isoformat(), 'source_commit':source_commit,
            'working_tree_dirty_at_run':dirty, 'source_sha256':source_hashes, 'config':config,
            'artifact_sha256':{str(p.relative_to(ROOT)):sha(p) for p in out.rglob('*') if p.is_file()},
            'software':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__},
            'hardware':platform.platform(),'elapsed_seconds':time.monotonic()-start,
            'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'limitations':['Four patients only; no publication-level statistical claim',
                          f'{len(genes)} shared biological targets after control exclusion and identifier harmonization of 541 assay features; 1000/3000-gene panels unavailable',
                          'Generic ImageNet ResNet18 reference, not pathology foundation embeddings',
                          'No validated cell/nucleus labels; RGB and texture summaries only',
                          'Simulated square/circular masks; tissue-area coverage unavailable',
                          'Thumbnail alignment audit cannot establish cell-level registration accuracy']}
    lines = ['# C0 development checkpoint', '', f"Status: **{spec['status']}**",'',
             f"Source commit: `{source_commit}`. Exact source hashes are in C0_SPEC_LOCK.json.", '',
             f"Prepared {sum(len(s.counts) for s in slides.values()):,} paired locations across four HEST IDC patients, {len(genes)} genes, and {len(mask_records)} frozen observation masks.", '',
             '## Verified', '',
             '- Unique barcode joins and patch/expression coordinate agreement within 1.01 full-resolution pixels.',
             '- Patient-disjoint train/validation/test assignments before target preprocessing.',
             '- Gene identifiers intersected and reordered by name; held-out expression values never used for panel ranking.',
             '- Training-only log-count normalization and gene ranking; raw counts retained separately.',
             '- Every frozen mask passed hidden-label NaN corruption invariance.',
             '- Raw assets and feature/target caches SHA-256 checked; feature barcode order verified.',
             '', '## Limitations', '']
    lines += ['- '+item for item in spec['limitations']]
    lines += ['', '## Next checkpoint', '',
              'C1 is not executed. Review C0 cohort/split suitability and encoder choice before model development, as specified by the supplied execution plan.',
              'The one-gene study is the top training log-variance gene. The 50/250 panels are available; the 1000-gene request is explicitly unsupported by this assay.',
              'No expression prediction, performance, biological discovery, or publication claims are made.', '',
              '## Reproduction', '', 'See README.md. All masks are in C0_MASKS/ and their geometry/checksums in C0_MASK_INDEX.json.',
              'C0_SPEC_LOCK.json refuses silent replacement. Later revisions must use a separately versioned checkpoint.']
    (out/'C0_REPORT.md').write_text('\n'.join(lines)+'\n')
    spec['artifact_sha256'] = {str(p.relative_to(ROOT)):sha(p) for p in out.rglob('*') if p.is_file()}
    write_json(out/'C0_SPEC_LOCK.json',spec)
    print(spec['status'], flush=True)
    if not passed:
        raise SystemExit(2)

if __name__ == '__main__':
    main()
