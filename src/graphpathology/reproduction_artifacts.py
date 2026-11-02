"""Integrity plus actual real-patient replay of recovered states, without refitting."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits

from .assets import sha
from .core import mlp_model, predict
from .reproduction_primary import checked_namespace, write_json
from .reproduction_replay import compare_tables, PANEL_METRICS


def matched_prediction(root, patient, item):
    from .backbone_train import load, design
    data = load(root, item['backbone'], patient)
    scaler = root / item['feature_scaler']; cp = root / item['checkpoint']
    if sha(scaler) != item['feature_scaler_sha256'] or sha(cp) != item['checkpoint_sha256']:
        raise ValueError('Recovered matched identity mismatch')
    with np.load(scaler, allow_pickle=False) as f:
        x = ((data['embeddings'] - f['mean']) / f['scale']).astype(np.float32)
    features, _ = design({patient:data}, {patient:x}, item['model'], item['control'], item['seed'])
    state = torch.load(cp, map_location='cpu', weights_only=True)
    if state['seed'] != item['seed'] or state['config'] != item['config']:
        raise ValueError('Matched saved configuration/seed mismatch')
    net = mlp_model(features[patient].shape[1], 50, state['config'])
    net.load_state_dict(state['state_dict'], strict=True); net.eval()
    image = predict(net, features[patient], 'B2').astype(float)
    if patient != 'TENX99':
        cal = json.loads((root/'checkpoints/legacy/C9/C9A/T1_SCALER.json').read_text())
        np.testing.assert_array_equal(data['genes'], cal['genes'])
        image = np.array(cal['destination']['center']) + np.array(cal['destination']['scale'])*(image-np.array(cal['source']['center']))/np.array(cal['source']['scale'])
    return data, image



def matched_reference(root, patient, run_id):
    """Reused states retain their old export path; unselected trials may lack exports."""
    base = root / 'experiments/journal_revision'
    options = [base / 'backbone_graph_comparison' / stage / patient / run_id / 'PANEL.csv.gz'
               for stage in ['primary_evaluation', 'evaluation']]
    options.append(base / 'transformer_egnn/evaluation' / patient / ('frozen_' + run_id) / 'PANEL.csv.gz')
    return next((p for p in options if p.is_file()), None)



def reference_threads(reference, default=4):
    """Read immutable runtime metadata; never choose a runtime from metric agreement."""
    if reference is None:
        return default
    metadata = reference.parent / 'COMPLETE.json'
    if not metadata.exists():
        return default
    return int(json.loads(metadata.read_text()).get('threads', default))


def verify_recovered(root, output, patient='BASS_P19'):
    from .extension_evaluate import model_inventory, load_data, prediction, masks, score_set
    output = checked_namespace(root, output, 'reproduction/runs')
    if output.exists(): raise FileExistsError(output)
    output.mkdir(parents=True); started = time.monotonic()
    data = load_data(root, patient)
    if len(data['targets']) > 2048: raise ValueError('Choose a genuinely measured small section for bounded verification')
    mask = next(m for m in masks(root,data,patient) if m['design']=='historical_V1' and m['budget']==.2 and m['mask_seed']==0)
    records = []; unique = set(); rows = []
    for item in model_inventory(root):
        cp = root / item['checkpoint']
        if sha(cp) != item['checkpoint_sha256']: raise ValueError('Recovered checkpoint identity mismatch')
        state = torch.load(cp,map_location='cpu',weights_only=True)
        if state['seed'] != item['seed']: raise ValueError('Saved seed mismatch')
        if item['origin']=='new':
            fit = json.loads((root/item['fit_record']).read_text())
            if state['epoch'] != fit['selected_epoch'] or state['variant'] != item['variant']:
                raise ValueError('Advanced selected-state link mismatch')
        elif state['config'] != item['old_config']: raise ValueError('Minibatch architecture mismatch')
        expected_path = root/'experiments/journal_revision/transformer_egnn/evaluation'/patient/item['id']/'PANEL.csv.gz'
        runtime = reference_threads(expected_path)
        torch.set_num_threads(runtime)
        image = prediction(root,data,patient,item)
        if image.shape!=(len(data['targets']),50) or not np.isfinite(image).all(): raise ValueError('Nonfinite recovered forward')
        panels,_,_,checks = score_set(data,image,mask,item,patient)
        expected_path = root/'experiments/journal_revision/transformer_egnn/evaluation'/patient/item['id']/'PANEL.csv.gz'
        expected = pd.read_csv(expected_path)
        expected = expected[expected.design.eq('historical_V1')&expected.budget.eq(.2)&expected.mask_seed.eq(0)]
        regression = compare_tables(pd.DataFrame(panels), expected, ['route'], PANEL_METRICS)
        records.append({'study':'advanced_48' if item['origin']=='new' else 'minibatch_30','identity':item['id'],
                        'checkpoint':item['checkpoint'],'sha256':sha(cp),'safe_load':True,'real_forward_rows':len(image),
                        'patient':patient,'torch_threads':runtime,'regression':regression,'hidden_corruption_checks':checks})
        unique.add(item['checkpoint']); rows.extend(panels)
        write_json(output/'ARTIFACT_CHECKS.json',records)
    candidates = json.loads((root/'experiments/journal_revision/backbone_graph_comparison/ALL_CANDIDATE_FIT_MANIFEST.json').read_text())
    selected = {r['run_id'] for r in json.loads((root/'experiments/journal_revision/backbone_graph_comparison/SELECTED_FIT_MANIFEST.json').read_text())}
    for item in candidates:
        expected_path = matched_reference(root,patient,item['run_id'])
        runtime = reference_threads(expected_path)
        torch.set_num_threads(runtime)
        data, image = matched_prediction(root,patient,item)
        if not np.isfinite(image).all(): raise ValueError('Nonfinite matched forward')
        meta = {'id':item['run_id'],'origin':'matched_'+item['backbone'],'variant':item['model']+'_'+item['control'],'seed':item['seed']}
        panels,_,_,checks = score_set(data,image,mask,meta,patient)
        expected_path = matched_reference(root,patient,item['run_id'])
        if expected_path is None:
            regression = {'status':'PARTIAL','scope':'Unselected candidate has no historical external export; safe real forward only.'}
        else:
            expected = pd.read_csv(expected_path)
            expected = expected[expected.design.eq('historical_V1')&expected.budget.eq(.2)&expected.mask_seed.eq(0)]
            regression = compare_tables(pd.DataFrame(panels),expected,['route'],PANEL_METRICS)
        records.append({'study':'matched_candidates_72','identity':item['run_id'],'checkpoint':item['checkpoint'],
                        'sha256':sha(root/item['checkpoint']),'safe_load':True,'real_forward_rows':len(image),'patient':patient,
                        'selected':item['run_id'] in selected,'torch_threads':runtime,'regression':regression,'hidden_corruption_checks':checks})
        unique.add(item['checkpoint']); rows.extend(panels)
        write_json(output/'ARTIFACT_CHECKS.json',records)
    # Optimizer/resume states are additional artifact files, never counted as fitted models.
    resumes=[]
    for path in sorted((root/'checkpoints/journal_revision/transformer_egnn').glob('*/resume.pt')):
        payload=torch.load(path,map_location='cpu',weights_only=True)
        if not {'model','optimizer','epoch'}<=payload.keys(): raise ValueError('Incomplete optimizer/resume state')
        resumes.append({'path':str(path.relative_to(root)),'sha256':sha(path),'epoch':payload['epoch'],'safe_load':True})
    write_json(output/'ARTIFACT_CHECKS.json',records);write_json(output/'RESUME_CHECKS.json',resumes)
    pd.DataFrame(rows).to_csv(output/'REAL_SAMPLE_PANEL.csv.gz',index=False)
    summary={'status':'EXACT_REPLAY' if all(r['regression']['status']=='EXACT_REPLAY' for r in records) else 'PARTIAL',
             'references_checked':len(records),'unique_best_states':len(unique),'safe_resume_files':len(resumes),
             'real_patient':patient,'locations':len(data['targets']),'real_mask':'original V1 random seed0,20% anchors',
             'scope':'All recovered model references; one actual small patient/mask, three prediction routes. No whole-cohort accuracy claim.',
             'reference_replays':sum(r['regression']['status']=='EXACT_REPLAY' for r in records),
             'safe_forward_only':sum(r['regression']['status']=='PARTIAL' for r in records),
             'elapsed_seconds':time.monotonic()-started}
    write_json(output/'SUMMARY.json',summary);return summary


def replay_matched(root, output):
    from .extension_evaluate import masks, score_set
    from .backbone_analysis import paired_effects
    from .reproduction_replay import patients
    output=checked_namespace(root,output,'reproduction/runs')
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True);started=time.monotonic();rows=[];regressions=[]
    items=json.loads((root/'experiments/journal_revision/backbone_graph_comparison/SELECTED_FIT_MANIFEST.json').read_text())
    items=[r for r in items if r['model'] in ['local','GCN'] and r['control'] in ['local','spatial']]
    if len(items)!=12:raise ValueError('Expected validation-selected local/GCN ×two backbones×three seeds')
    for cohort,patient in patients(root):
        if cohort!='bassiouni':continue
        for item in items:
            expected_path=matched_reference(root,patient,item['run_id'])
            runtime=reference_threads(expected_path)
            torch.set_num_threads(runtime)
            data,image=matched_prediction(root,patient,item)
            ms=[m for m in masks(root,data,patient) if m['design']=='historical_V1' and m['budget']==.2]
            if len(ms)!=3:raise ValueError('Expected three fixed neural mask seeds')
            panels=[]
            for mask in ms:
                meta={'id':item['run_id'],'origin':'matched_'+item['backbone'],'variant':item['model']+'_'+item['control'],'seed':item['seed']}
                panel,_,_,_=score_set(data,image,mask,meta,patient)
                panels.extend([r for r in panel if r['route']=='exploratory_residual_anchor'])
            expected_path=matched_reference(root,patient,item['run_id'])
            if expected_path is None:raise FileNotFoundError('Selected model must have historical external exports')
            expected=pd.read_csv(expected_path)
            expected=expected[expected.design.eq('historical_V1')&expected.budget.eq(.2)&expected.route.eq('exploratory_residual_anchor')]
            actual=pd.DataFrame(panels)
            regression=compare_tables(actual,expected,['mask_seed'],PANEL_METRICS)
            regressions.append({'patient':patient,'run_id':item['run_id'],**regression})
            actual=actual.assign(backbone=item['backbone'],model=item['model'])
            rows.append(actual)
        print('MATCHED REPLAYED',patient,flush=True)
    panel=pd.concat(rows,ignore_index=True);panel.to_csv(output/'PANEL.csv.gz',index=False)
    # Three equal masks within seed, then three equal training seeds within patient.
    keys=['backbone','patient','model'];patient=panel.groupby(keys+['training_seed'])[PANEL_METRICS].mean().reset_index().groupby(keys)[PANEL_METRICS].mean().reset_index()
    effects=[]
    for backbone,g in patient.groupby('backbone'):
        a=g[g.model.eq('GCN')].set_index('patient');b=g[g.model.eq('local')].set_index('patient')
        for pid in a.index:
            effects.append({'backbone':backbone,'patient':pid,'delta_pearson':a.loc[pid,'median_pearson']-b.loc[pid,'median_pearson'],
                            'relative_mse':a.loc[pid,'mean_mse']/b.loc[pid,'mean_mse']-1})
    effects=pd.DataFrame(effects);effects.to_csv(output/'PAIRED_PATIENT_EFFECTS.csv',index=False)
    effects.groupby('backbone')[['delta_pearson','relative_mse']].median().to_csv(output/'COHORT_EFFECTS.csv')
    write_json(output/'REGRESSION.json',regressions)
    summary={'status':'EXACT_REPLAY' if all(r['status']=='EXACT_REPLAY' for r in regressions) else 'PARTIAL','patients':22,
             'selected_models':12,'panel_rows':len(panel),'elapsed_seconds':time.monotonic()-started,
             'scope':'Actual recovered-weight forward predictions, original three neural masks,20% distributed residual GCN-versus-local, both backbones.'}
    write_json(output/'SUMMARY.json',summary);return summary


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('action',choices=['verify','matched']);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
    with threadpool_limits(4):result=verify_recovered(a.root.resolve(),a.output.resolve()) if a.action=='verify' else replay_matched(a.root.resolve(),a.output.resolve())
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
