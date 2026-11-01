import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from support import *
from policy import choose,region

def test_sampling_ignores_expression_and_matches_budgets():
 xy=np.random.default_rng(5).normal(size=(120,2));forbidden=region(xy,'V2',1);eligible=np.flatnonzero(~forbidden)
 for strategy in ['random','spatial_FPS']:
  order=choose({'xy':xy},eligible,48,strategy,3)
  for values in [np.full((120,280),np.nan),np.random.default_rng(1).normal(size=(120,280))]:assert np.array_equal(order,choose({'xy':xy,'expression':values},eligible,48,strategy,3))
  assert len(np.unique(order))==48 and not forbidden[order].any()
  for b in [.1,.2,.4]:assert len(order[:count(120,b)])==count(120,b)

def test_frozen_models_corruption_and_gene_restriction():
 for panel in panels():
  s=scene('A1',panel);y=np.random.default_rng(8).normal(size=s['image'].shape).astype('float32');obs=np.arange(len(y))%5==0;h=np.flatnonzero(~obs)
  for model in ['E0','E1','E2','E3']:
   p=predict(s,c5.observed_payload(y,obs,s['xy']),h,panel,'V1',model)
   for val in [np.nan,12345.]:
    alt=y.copy();alt[~obs]=val;q=predict(s,c5.observed_payload(alt,obs,s['xy']),h,panel,'V1',model);assert np.array_equal(p,q)
   assert p.shape==(len(h),len(panels()[panel])) and np.isfinite(p).all()

def test_source_coordinate_identity():
 for sample in SAMPLES:
  with np.load(CACHE/(sample+'_metadata.npz')) as m,np.load(CACHE/(sample+'_features.npz')) as f:
   assert np.array_equal(m['barcodes'],f['barcodes']) and np.array_equal(m['xy_um'],f['xy_um']);assert len(np.unique(m['barcodes']))==len(m['barcodes']);assert len(m['barcodes'])>=100
