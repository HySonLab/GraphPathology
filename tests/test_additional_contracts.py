from pathlib import Path
import json
import numpy as np
import pytest
from graphpathology import policy,core,external
ROOT=Path(__file__).resolve().parents[1]
LEGACY=ROOT/'checkpoints/legacy'

def test_original_search_budgets_and_topology_schedules():
    c=json.loads((LEGACY/'C4/C4A/configs/protocol.json').read_text())
    assert len(c['encoders'])==2 and len(c['mlp_candidates'])==4 and len(c['ridge_alphas'])==5 and c['training_seeds']==[0,1,2]
    c=json.loads((LEGACY/'C3/C3A/configs/protocol.json').read_text())
    assert c['training_seeds']==[0,1,2] and c['graph_randomization_seeds']==[0,1,2,3,4]
    assert c['conditional_lopo_new_controls']==['A1','A6']

def test_calibration_folds_and_selection_no_outer_access():
    from graphpathology.selection import select,require_training_patients
    c=json.loads((LEGACY/'C8/configs/protocol.json').read_text())
    assert sorted(sum(c['folds'].values(),[]))==[x+'1' for x in 'ABCDEFGH']
    for outer in c['folds'].values():
        train=sorted(set(x+'1' for x in 'ABCDEFGH')-set(outer))
        assert len(require_training_patients(train,outer))==6
        with pytest.raises(ValueError):require_training_patients(outer,outer)
    rows=[]
    for s in ['A1','B1']:
        rows.extend([{'candidate':'T0','sample_id':s,'regime':'V1','budget':.2,'median_pearson':.2,'common_mse':2.},{'candidate':'T1','sample_id':s,'regime':'V1','budget':.2,'median_pearson':.204,'common_mse':3.}])
    assert select(rows,'candidate')[0]=='T0'

def test_external_sampling_uses_no_expression():
    xy=np.random.default_rng(5).normal(size=(120,2));f=policy.region(xy,'V2',1);eligible=np.flatnonzero(~f)
    for strategy in ['random','spatial_FPS']:
        a=policy.choose({'xy':xy},eligible,48,strategy,3)
        for v in [np.full((120,50),np.nan),np.random.default_rng(4).normal(size=(120,50))]:np.testing.assert_array_equal(a,policy.choose({'xy':xy,'expression':v},eligible,48,strategy,3))
        assert len(set(a))==48 and not f[a].any()

def test_missing_original_heads_error_is_specific(tmp_path):
    from graphpathology.workflow import infer
    cfg={'max_locations':20,'models':[f'C4/C4A/models/C4A_primary_resnet18_reference_primary_S2-50_D1_h256d0.1_s{i}.pt' for i in range(3)]}
    p=tmp_path/'config.json';p.write_text(json.dumps(cfg))
    with pytest.raises(FileNotFoundError,match='4aaa11ab8998627d'):infer(ROOT,p,True)

def test_migration_small_records_and_prefix_scope(tmp_path):
    import csv,hashlib
    from graphpathology.assets import resolve,verify
    (tmp_path/'checkpoints/legacy/results').mkdir(parents=True)
    old='results/registry.csv';new='checkpoints/legacy/'+old;data=b'original prefix\nappended row\n'
    (tmp_path/new).write_bytes(data)
    for name,rows in [('path_map.csv',[dict(old_path=old,new_path=new,sha256=hashlib.sha256(data).hexdigest())]),('manifest.csv',[dict(original_path=old,path=new)])]:
        with (tmp_path/'checkpoints'/name).open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    row=dict(path=old,expected_sha256=hashlib.sha256(data[:16]).hexdigest(),scope='prefix',prefix_bytes=16,declared_in='fixture',starting_status='match')
    with (tmp_path/'checkpoints/historical_expectations.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=row);w.writeheader();w.writerow(row)
    assert resolve(tmp_path,old)==tmp_path/new
    assert verify(tmp_path,True)['historical_counts']=={'match':1}
    (tmp_path/new).write_bytes(b'changed')
    assert verify(tmp_path)['migration']['changed']==1

def test_c8_gate_families_hidden_invariance_fixture():
    from graphpathology.gating import features,predict_B
    rng=np.random.default_rng(12);xy=rng.normal(size=(40,2))*300;y=rng.normal(size=(40,40));im=rng.normal(size=(40,40));obs=np.arange(0,40,4);h=np.setdiff1d(np.arange(40),obs)
    obj={'source':external.moments([y[:10]]),'destination':external.moments([y[10:20]*2+1])}
    p,dist=external.predict(xy,im,y,obs,h,'V1',obj)
    def payload(target):
        predictions,_=external.predict(xy,im,target,obs,h,'V1',obj)
        center=np.array(obj['destination']['center']);scale=np.array(obj['destination']['scale'])
        return {'image':(im[h]-obj['source']['center'])/obj['source']['scale'],'st':(predictions['M1']-center)/scale,'features':features(xy,obs,h,target,obj['destination']),'distance':dist}
    before=payload(y);corrupt=y.copy();corrupt[h]=np.nan;after=payload(corrupt)
    for model in [{'family':'F0'},{'family':'F1','lambda':.2},{'family':'F2','weights':[0,.2,.5]},{'family':'F3','feature_mean':[0.]*4,'feature_scale':[1.]*4,'coef':[-1.,.1,.2,.3,.4]}]:np.testing.assert_array_equal(predict_B(before,model),predict_B(after,model))

def test_c10_frozen_c9_scaler_and_no_new_fitting():
    cfg=json.loads((LEGACY/'C10/configs/protocol.json').read_text());scaler=json.loads((LEGACY/'C9/C9A/T1_SCALER.json').read_text())
    assert cfg['calibration_patients']==[]
    assert not set(cfg['inputs'])&set(scaler['calibration_patients'])
    assert cfg['genes']==scaler['genes']

def test_table_s8_schema_matches_production_parameter_counts():
    from graphpathology.graph import parameter_count
    rows=json.loads((ROOT/'manuscript/history/tables/TableS8.json').read_text())['rows']
    assert len(rows)==9 and len({r['predictor'] for r in rows})==9
    for row in rows:
        expected=parameter_count(row['input'],row['hidden'],row['output']) if row['hidden'] else (row['input']+1)*row['output']
        assert row['parameters']==expected
