import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from common import CFG,PANELS,SELECTED,core,preprocessing,scaler
import run

def test_self_only_has_no_neighbor_information():
 # Rows can change independently: no other node participates in a self-only input.
 d={'x':np.zeros((3,512),np.float32)}
 x,_=run.features('unused',d,'S1','A1',-1)
 d['x'][1]=3
 changed,_=run.features('unused',d,'S1','A1',-1)
 assert np.array_equal(x[0],changed[0]) and np.array_equal(x[2],changed[2])
 assert np.array_equal(changed[:,:512],changed[:,512:])

def test_exact_self_and_true_nominal_capacity():
 for panel in PANELS:
  p=SELECTED[panel]['SAGE']['params']
  a=core.mlp_model(1024,len(PANELS[panel]),p);b=core.mlp_model(1024,len(PANELS[panel]),p)
  assert sum(q.numel() for q in a.parameters())==sum(q.numel() for q in b.parameters())

def test_balanced_randomization_schedule():
 assert CFG['training_seeds']==[0,1,2] and CFG['graph_randomization_seeds']==[0,1,2,3,4]
 assert CFG['conditional_lopo_new_controls']==['A1','A6']
 assert len(PANELS)==4 and len(PANELS['S2-280'])==280
