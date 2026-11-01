import numpy as np
import pytest
from graphpathology.data import validate_patient_splits
from graphpathology.preprocessing import CountTransform, select_panels
from graphpathology.masks import make_mask, anchor_distances, stable_seed


def test_hidden_labels_never_reach_transform():
    counts = np.arange(60).reshape(12,5)
    transform = CountTransform().fit(counts[:4])
    mask = np.arange(12)%3==0
    original,_ = transform.observed_inputs(counts,mask)
    corrupt = counts.astype(float); corrupt[~mask]=np.nan
    altered,_ = transform.observed_inputs(corrupt,mask)
    np.testing.assert_array_equal(original,altered)
    assert np.all(original[~mask]==0)
    assert np.all(transform.observed_inputs(corrupt,np.zeros(12,dtype=bool))[0]==0)


def test_heldout_targets_cannot_change_fitted_parameters_or_panel():
    train = np.array([[0,1,4],[1,1,3],[8,1,1]])
    transform=CountTransform().fit(train)
    expected=transform.mean.copy()
    transform.transform(np.full((10,3),1000000))
    np.testing.assert_array_equal(transform.mean,expected)
    panel=select_panels(train,['a','b','c'],(1,1000))
    assert panel['1']['genes']==['a']
    assert panel['1000']['status']=='unavailable_panel_too_small'


def test_patient_leakage_rejected():
    with pytest.raises(ValueError,match='leakage'):
        validate_patient_splits({'train':['a'],'validation':['b'],'test':['c']}, {'a':'p1','b':'p1','c':'p2'})
    validate_patient_splits({'train':['a'],'validation':['b'],'test':['c']}, {'a':'p1','b':'p2','c':'p3'})


@pytest.mark.parametrize('regime',['random','block','capture_window'])
@pytest.mark.parametrize('q',[.05,.1,.2,.4,.6,.8])
def test_frozen_mask_geometry_and_distance(regime,q):
    xy=np.stack(np.meshgrid(np.arange(30)*100,np.arange(30)*100),axis=-1).reshape(-1,2)
    seed=stable_seed('slide',regime,q,0)
    mask,meta=make_mask(xy,regime,q,seed)
    np.testing.assert_array_equal(mask,make_mask(xy,regime,q,seed)[0])
    assert mask.any() and (~mask).any()
    if regime!='random':
        delta=xy-np.array(meta['center_um'])
        distance=np.max(abs(delta),axis=1) if regime=='block' else np.linalg.norm(delta,axis=1)
        inside=distance<=meta['radius_um']
        np.testing.assert_array_equal(mask,~inside if regime=='block' else inside)
    else:
        assert mask.sum()==round(len(xy)*q)
    d=anchor_distances(xy,mask)
    assert np.all(d[mask]==0) and np.all(d[~mask]>0)


def test_observed_zeros_have_distinct_missing_flag():
    transform=CountTransform().fit(np.array([[0,0],[2,3]]))
    _,flag=transform.observed_inputs(np.zeros((2,2)),np.array([True,False]))
    assert flag.tolist()==[True,False]


def test_invalid_shapes_and_coverage():
    with pytest.raises(ValueError): make_mask(np.zeros((10,2)),'random',1,0)
    with pytest.raises(ValueError): CountTransform().fit(np.array([[-1]]))
    with pytest.raises(ValueError): CountTransform().fit(np.ones((2,2))).observed_inputs(np.ones((2,2)),np.ones(2))
