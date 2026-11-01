import importlib.util
from pathlib import Path
import numpy as np
import pytest
from graphpathology import core as c

def test_constant_baseline_correlations_undefined():
 y=np.arange(30).reshape(10,3);pred=np.full((10,3),.1)
 m=c.gene_metrics(y,pred)
 assert np.isnan(m['pearson']).all() and np.isnan(m['spearman']).all()
 assert c.coexpression(y,pred)[0]['n_valid_pairs']==0

def test_metrics_known_answer_and_coexpression():
 rng=np.random.default_rng(4);y=rng.normal(size=(30,4));p=2*y+3
 np.testing.assert_allclose(c.gene_metrics(y,p)['pearson'],1)
 assert c.coexpression(y,p)[0]['coexpression_pearson']==pytest.approx(1)
 np.testing.assert_allclose(c.gene_metrics(y,y)['mse'],0)

def test_scaler_training_only_roundtrip_and_heldout_invariance():
 a=np.arange(40).reshape(10,4);s=c.Scaler().fit(a);before=s.mean.copy()
 s.transform(np.full((100,4),1e9));np.testing.assert_array_equal(before,s.mean)
 np.testing.assert_allclose(s.inverse(s.transform(a)),a,atol=1e-5)

def test_lopo_secondary_scaling_cancels_primary_affine_transform():
 rng=np.random.default_rng(0);a=rng.normal(size=(30,5));shift=np.arange(5);scale=np.arange(5)+2
 primary=(a-shift)/scale
 local1=c.Scaler().fit(a[:20]);local2=c.Scaler().fit(primary[:20])
 np.testing.assert_allclose(local1.transform(a),local2.transform(primary),atol=1e-6)

def test_negative_controls_preserve_validation_and_gene_vectors():
 x=np.arange(60).reshape(20,3);y=np.arange(80).reshape(20,4)
 xp,yp,index=c.permute_training(x,y,'NC1',2)
 np.testing.assert_array_equal(xp,x);np.testing.assert_array_equal(yp,y[index])
 xp,yp,index=c.permute_training(x,y,'NC2',2)
 np.testing.assert_array_equal(xp,x[index]);np.testing.assert_array_equal(yp,y)

def test_nonfinite_prediction_rejected():
 with pytest.raises(ValueError): c.gene_metrics(np.ones((3,2)),np.full((3,2),np.nan))


def test_mlp_seed_reproducibility_and_fixed_epoch_without_validation():
 c.torch.set_num_threads(1)
 rng=np.random.default_rng(6);x=rng.normal(size=(20,3)).astype(np.float32);y=x[:,:2]
 cfg={'hidden':8,'dropout':.1};training={'learning_rate':.001,'weight_decay':.0001,'batch_size':8,'max_epochs':3,'patience':2,'min_delta':.0001}
 a,info=c.train_mlp(x,y,None,None,cfg,training,2,fixed_epochs=3)
 b,_=c.train_mlp(x,y,None,None,cfg,training,2,fixed_epochs=3)
 np.testing.assert_array_equal(c.predict(a,x,'B2'),c.predict(b,x,'B2'))
 assert info['best_epoch']==3


def test_ridge_single_gene_keeps_matrix_contract():
 rng=np.random.default_rng(9);x=rng.normal(size=(30,4));y=x[:,:1]
 model=c.Ridge(alpha=1).fit(x,y)
 p=c.predict(model,x,'B1')
 assert p.shape==(30,1)
 assert c.gene_metrics(y,p)['pearson'][0]>.99
