import h5py
import numpy as np
import pytest
from graphpathology.data import read_slide


def fixture(root):
    (root/'IDC/adata').mkdir(parents=True)
    (root/'IDC/patches').mkdir()
    with h5py.File(root/'IDC/adata/demo.h5ad','w') as f:
        f['X']=np.array([[1,2],[3,4],[5,6]])
        f['obs/_index']=np.array([b'a',b'b',b'c'])
        f['var/_index']=np.array([b'g1',b'g2'])
        f['obsm/spatial']=np.array([[10.3,20.1],[30.2,40.5],[50.6,60.7]])
        f['uns/spatial/ST/images/downscaled_fullres']=np.zeros((100,100,3),dtype=np.uint8)
        f['uns/spatial/ST/scalefactors/tissue_downscaled_fullres_scalef']=1.
    with h5py.File(root/'IDC/patches/demo.h5','w') as f:
        f['barcode']=np.array([[b'c'],[b'a']])
        f['coords']=np.array([[60,50],[20,10]])
        f['img']=np.zeros((2,224,224,3),dtype=np.uint8)


def test_barcode_join_handles_reordered_subset(tmp_path):
    fixture(tmp_path)
    slide=read_slide(tmp_path,'demo',.5)
    np.testing.assert_array_equal(slide.counts,[[5,6],[1,2]])
    np.testing.assert_allclose(slide.xy_um,[[25.3,30.35],[5.15,10.05]])
    assert slide.excluded_expression_locations==1


def test_wrong_coordinate_orientation_rejected(tmp_path):
    fixture(tmp_path)
    with h5py.File(tmp_path/'IDC/patches/demo.h5','r+') as f:
        f['coords'][:]=f['coords'][:][:,::-1]
    with pytest.raises(ValueError,match='coordinate mismatch'):
        read_slide(tmp_path,'demo',.5)


def test_duplicate_barcodes_rejected(tmp_path):
    fixture(tmp_path)
    with h5py.File(tmp_path/'IDC/patches/demo.h5','r+') as f:
        f['barcode'][:]=np.array([[b'a'],[b'a']])
    with pytest.raises(ValueError,match='Duplicate'):
        read_slide(tmp_path,'demo',.5)


def test_assay_controls_are_not_prediction_targets(tmp_path):
    fixture(tmp_path)
    with h5py.File(tmp_path/'IDC/adata/demo.h5ad','r+') as f:
        del f['var/_index']
        f['var/_index']=np.array([b'GENE',b'NegControlProbe_001'])
    slide=read_slide(tmp_path,'demo',.5)
    assert slide.genes.tolist()==['GENE']
    assert slide.excluded_control_features==['NegControlProbe_001']
    assert slide.counts.shape==(2,1)


def test_harmonization_joins_gene_names_not_column_position():
    from types import SimpleNamespace
    from graphpathology.data import harmonize_genes
    slides={'a':SimpleNamespace(genes=np.array(['g2','g1','g3']),counts=np.array([[20,10,30]])),
            'b':SimpleNamespace(genes=np.array(['g1','g2']),counts=np.array([[100,200]]))}
    record=harmonize_genes(slides)
    np.testing.assert_array_equal(slides['a'].counts,[[10,20]])
    np.testing.assert_array_equal(slides['b'].counts,[[100,200]])
    assert record['a']['nonshared_biological_features']==['g3']
