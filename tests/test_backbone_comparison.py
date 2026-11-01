import numpy as np
import pandas as pd
import pytest
from graphpathology.backbone_analysis import paired_effects,common_valid_runs
from graphpathology.journal_learned import hybrids
from graphpathology.graph import build_graph,shuffle_graph


def test_relative_error_is_paired_before_cohort_summary():
    base={'cohort':['Bassiouni']*2,'patient':['p1','p2'],'route':['residual']*2,'design':['V1']*2,'condition':['distributed']*2,'budget':[.2]*2,'median_pearson':[.4,.5],'mean_mae':[1.,1.]}
    ref=pd.DataFrame(base|{'mean_mse':[1.,10.]});candidate=pd.DataFrame(base|{'mean_mse':[2.,9.]})
    result=paired_effects(candidate,ref)
    assert np.median(result.relative_mse)==pytest.approx(.45)
    assert np.median(candidate.mean_mse)/np.median(ref.mean_mse)-1==0


def test_common_valid_gene_sensitivity_does_not_substitute_primary_metric():
    base={'patient':['p']*3,'design':['V1']*3,'condition':['distributed']*3,'budget':[.2]*3,'mask_seed':[0]*3,'geometry':[0]*3,'training_seed':[0]*3,'mask_sha256':['same']*3,'gene':['a','b','c']}
    candidate=pd.DataFrame(base|{'pearson':[.9,.1,np.nan]});reference=pd.DataFrame(base|{'pearson':[.8,np.nan,.2]})
    result=common_valid_runs(candidate,reference)
    assert result.common_valid_genes.iloc[0]==1
    assert result.common_delta_pearson.iloc[0]==pytest.approx(.1)
    assert candidate.pearson.median()-reference.pearson.median()==0
    reference['mask_sha256']='different'
    with pytest.raises(ValueError,match='mismatch'):common_valid_runs(candidate,reference)


def test_hybrid_corrupts_unscored_unobserved_labels_too():
    rng=np.random.default_rng(3);xy=rng.normal(size=(64,2));y=rng.normal(size=(64,5));image=rng.normal(size=(64,5));observed=np.arange(16);hidden=np.arange(16,32)
    original=hybrids(xy,y,image,observed,hidden)
    for mode in ['nan','random']:
        corrupt=y.copy();corrupt[16:]=np.nan if mode=='nan' else rng.normal(size=(48,5))
        result=hybrids(xy,corrupt,image,observed,hidden)
        for route in original:np.testing.assert_array_equal(result[route],original[route])
    alternative=hybrids(xy,y,image,np.arange(32,48),hidden)
    np.testing.assert_array_equal(alternative['image'],original['image'])


def test_shuffled_control_changes_spatial_alignment():
    xy=np.random.default_rng(4).uniform(0,1000,(200,2));a=build_graph(xy,{'kind':'delaunay','max_edge_um':250});b,p=shuffle_graph(a,'section',0)
    assert 1-a.multiply(b).nnz/a.nnz>.5
    np.testing.assert_array_equal(np.sort(p),np.arange(200))
    assert np.allclose(np.sort(np.asarray(a.sum(1)).ravel()),np.sort(np.asarray(b.sum(1)).ravel()))
