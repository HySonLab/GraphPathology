import numpy as np
import pytest
from scipy.spatial import cKDTree
from graphpathology.journal_methods import candidates, predict_anchors, identical_target_design, patient_bootstrap


def test_anchor_only_all_unobserved_corruption_and_original_idw():
    rng = np.random.default_rng(12); xy = rng.normal(size=(50, 2)); y = rng.normal(size=(50, 3))
    obs = np.array([3, 2, 30, 18, 22]); query = np.array([8, 9, 10, 11])
    corrupt = y.copy(); corrupt[np.setdiff1d(np.arange(50), obs)] = np.nan
    for config in candidates():
        p, _ = predict_anchors(xy, y, obs, query, config)
        q, _ = predict_anchors(xy, corrupt, obs, query, config)
        np.testing.assert_array_equal(p, q)
    d, j = cKDTree(xy[obs]).query(xy[query], k=5); w = 1 / np.maximum(np.round(d, 6), 1e-6) ** 2; w /= w.sum(1, keepdims=True)
    expected = np.einsum('nk,nkg->ng', w, y[obs][j], optimize=True)
    p, _ = predict_anchors(xy, y, obs, query, {'kind': 'idw', 'k': 16, 'p': 2})
    np.testing.assert_array_equal(expected, p)
    with pytest.raises(ValueError, match='overlap'): predict_anchors(xy, y, obs, obs[:1], candidates()[0])


def test_fixed_target_design_equal_counts_nested_and_label_free():
    xy = np.random.default_rng(8).normal(size=(1000, 2))
    r, h, orders, _ = identical_target_design(xy, 'P1', 0)
    r1, h1, orders1, _ = identical_target_design(xy, 'P1', 1)
    np.testing.assert_array_equal(r, r1); np.testing.assert_array_equal(h, h1)
    assert len(h) == 80 and r.sum() == 400 and r[h].all()
    for name, order in orders.items():
        assert not np.intersect1d(order, h).size
        assert len(order[:100]) == len(order[:400]) / 4
        if name == 'region_exclusion': assert not r[order].any()
    assert not np.array_equal(orders['distributed'], orders1['distributed'])


def test_patient_bootstrap_pairs_and_undefined():
    a = patient_bootstrap([1, 2, 3, np.nan])
    b = patient_bootstrap([1, 2, 3])
    assert a == b and a['n_patients'] == 3 and a['median_effect'] == 2 and a['fraction_improved'] == 1
    assert a['ci_low'] == 1 and a['ci_high'] == 3
    assert patient_bootstrap([-1, 0, 1], direction='negative')['fraction_improved'] == 1 / 3
    assert patient_bootstrap([np.nan])['n_patients'] == 0


def test_constant_anchors_and_extreme_gaussian_are_finite():
    xy = np.array([[0, 0], [0, 0], [1e10, 1e10]], float); y = np.ones((3, 2))
    for c in candidates():
        p, _ = predict_anchors(xy, y, [0, 1], [2], c)
        np.testing.assert_allclose(p, 1, atol=1e-14)


def test_patient_pairing_rejects_silent_drops():
    import pandas as pd
    from graphpathology.journal_statistics import contrast
    a=pd.DataFrame({'patient':['P1'],'median_pearson':[.2],'mean_mse':[1.],'mean_mae':[.5]})
    b=a.copy();b['patient']='P2'
    with pytest.raises(ValueError,match='Unmatched'):contrast(a,b,['patient'],'comparison','test')


def test_learned_hybrids_never_read_unobserved_labels():
    from graphpathology.journal_learned import hybrids
    rng=np.random.default_rng(4);xy=rng.normal(size=(30,2));y=rng.normal(size=(30,3));image=rng.normal(size=(30,3))
    observed=np.array([1,4,8,12]);hidden=np.array([9,10,11]);bad=y.copy();bad[np.setdiff1d(np.arange(30),observed)]=np.nan
    a=hybrids(xy,y,image,observed,hidden);b=hybrids(xy,bad,image,observed,hidden)
    for route in a:np.testing.assert_array_equal(a[route],b[route])
