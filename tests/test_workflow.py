import json,sys,subprocess
from pathlib import Path
import numpy as np
import pytest
from scipy import sparse
from graphpathology import external,policy,graph,core
ROOT=Path(__file__).resolve().parents[1]

def test_graph_counts_and_known_operators():
    assert graph.parameter_count(512,256,50)==144178
    assert graph.parameter_count(512,128,50)==72114
    assert graph.parameter_count(1024,67,50)==72075
    a=sparse.csr_matrix([[0,1,1,0],[1,0,0,0],[1,0,0,0],[0,0,0,0]])
    np.testing.assert_allclose(graph.operator(a,'MEAN').toarray(),[[1/3,1/3,1/3,0],[.5,.5,0,0],[.5,0,.5,0],[0,0,0,1]])
    np.testing.assert_allclose(graph.operator(a,'GCN').toarray(),[[1/3,1/np.sqrt(6),1/np.sqrt(6),0],[1/np.sqrt(6),.5,0,0],[1/np.sqrt(6),0,.5,0],[0,0,0,1]])
    np.testing.assert_allclose(graph.operator(a,'SAGE').toarray(),[[0,.5,.5,0],[1,0,0,0],[1,0,0,0],[0,0,0,0]])
    assert not np.allclose(graph.operator(a,'GCN').sum(1),1)

def test_knn_union_tie_by_index():
    xy=np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]],float)
    a=graph.build_graph(xy,{'kind':'knn','k':1,'max_edge_um':1})
    assert set(zip(*sparse.triu(a).nonzero()))=={(0,1),(0,2),(0,3),(0,4)}

def test_nested_masks_and_no_hidden_leakage_including_unscored():
    xy=np.array([(x,y) for x in range(8) for y in range(8)],float)
    rng=np.random.default_rng(1);y=rng.normal(size=(64,3));im=rng.normal(size=(64,3));obj={'source':external.moments([y[:8]]),'destination':external.moments([y[8:16]])}
    for geometry in range(3):
        f=policy.region(xy,'V2',geometry);eligible=np.flatnonzero(~f);order=policy.choose({'xy':xy},eligible,policy.count(64,.4),'random',2)
        assert f.sum()==26 and not f[order].any()
        sets=[set(order[:policy.count(64,b)]) for b in [.1,.2,.4]];assert sets[0]<sets[1]<sets[2]
        obs=order[:13];hidden=np.flatnonzero(f);p,_=external.predict(xy,im,y,obs,hidden,'V2',obj)
        unobserved=np.setdiff1d(np.arange(64),obs)
        for bad in [np.nan,912.]:
            corrupt=y.copy();corrupt[unobserved]=bad;q,_=external.predict(xy,im,corrupt,obs,hidden,'V2',obj)
            for method in p:np.testing.assert_array_equal(p[method],q[method])

def test_manual_interpolation_and_calibration_direction():
    xy=np.array([[0,0],[3,0],[1,0]],float);y=np.array([[2.],[8.],[999.]])
    obj={'source':{'center':[1.],'scale':[2.]},'destination':{'center':[10.],'scale':[3.]}}
    p,_=external.predict(xy,np.full((3,1),5.),y,np.array([0,1]),np.array([2]),'V1',obj)
    assert p['M1'][0,0]==pytest.approx(3.2)
    assert p['M3'][0,0]==pytest.approx(.8*3.2+.2*16)
    p,_=external.predict(xy,np.full((3,1),5.),y,np.array([0,1]),np.array([2]),'V2',obj)
    assert p['M1'][0,0]==5
    xy[2]=xy[0];p,_=external.predict(xy,np.ones((3,1)),y,np.array([0]),np.array([2]),'V1',obj)
    assert p['M1'][0,0]==2 and np.isfinite(p['M3']).all()

def test_gate_population_variance_in_destination_units():
    from graphpathology.gating import features
    xy=np.array([[0.,0],[100,0],[200,0]])
    y=np.array([[2.,4],[6,12],[np.nan,np.nan]])
    f=features(xy,np.array([0,1]),np.array([2]),y,{'center':[0,0],'scale':[2,4]})
    assert f[0,3]==pytest.approx(np.log1p(1.))
    assert f[0,2]==pytest.approx(np.log1p(2))

def test_residual_messages_design_and_intercept():
    from graphpathology import residual as r
    xy=np.array([[0.,0],[100,0],[200,0]])
    a=sparse.csr_matrix([[0,1,0],[1,0,1],[0,1,0]])
    r.GRAPHS['manual']=graph.operator(a,'MEAN')
    scene={'xy':xy,'sample':'manual','image':np.ones((3,1)),'projection':np.ones((3,8),dtype='float32')}
    anchors=r.observed_payload(np.array([[6.],[np.nan],[np.nan]]),np.array([True,False,False]),xy)
    parts=r.feature_parts(scene,anchors,{'best_family':'B2a','best_param':1},.2)
    shared,gene,_,_=parts
    assert shared.shape==(3,29) and r.design(shared,gene,0).shape==(3,35)
    np.testing.assert_allclose(shared[:,-1],1/16)
    np.testing.assert_allclose(shared[:,25],[.5,1/3,0])
    np.testing.assert_allclose(shared[:,26],[5/12,5/18,1/6])
    x=np.column_stack([np.zeros(3),np.ones(3)])
    np.testing.assert_allclose(r.fit_residual([x],[np.full(3,7.)]),[[0,7]])

def test_ties_and_paired_aggregation():
    y=np.array([[0.],[0.],[1.],[2.]])
    assert core.gene_metrics(y,y)['spearman'][0]==pytest.approx(1)
    a=np.array([0.,10,11]);b=np.array([0.,1,10])
    assert np.median(a-b)==1 and np.median(a)-np.median(b)==9

@pytest.mark.parametrize('command',['doctor','verify-archive','resolve','smoke','publication','build-documents','infer'])
def test_help_from_outside_checkout(command,tmp_path):
    result=subprocess.run([sys.executable,'-m','graphpathology',command,'--help'],cwd=tmp_path,capture_output=True,text=True,timeout=15)
    assert result.returncode==0 and 'usage:' in result.stdout

def test_smoke_serialization_and_output_safety(tmp_path):
    from graphpathology.workflow import smoke
    from graphpathology.assets import new_output
    output=tmp_path/'path with spaces'
    assert smoke(ROOT,output)['training_steps']==2
    assert (output/'fixture_metrics.csv').exists()
    with pytest.raises(FileExistsError):smoke(ROOT,output)
    with pytest.raises(ValueError):new_output(ROOT,ROOT/'checkpoints/legacy/new')

def test_frozen_encoder_dropout_and_three_head_inverse_scaling():
    import torch
    from graphpathology.vision import resnet18,ensemble,patch_transform
    from PIL import Image
    model=resnet18();assert not model.training and not any(p.requires_grad for p in model.parameters())
    assert patch_transform(Image.new('RGB',(12,15))).shape==(3,224,224)
    sc=core.Scaler().fit(np.array([[0.,0],[2,2]]));ts=core.Scaler().fit(np.array([[10.],[14.]]))
    heads=[core.mlp_model(2,1,{'hidden':2,'dropout':.1}) for _ in range(3)]
    for i,h in enumerate(heads):
        for p in h.parameters():p.data.zero_()
        h[-1].bias.data.fill_(i)
    np.testing.assert_array_equal(ensemble(heads,np.zeros((3,2)),sc,ts),np.full((3,1),14.))
