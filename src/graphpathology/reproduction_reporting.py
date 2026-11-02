"""Assemble new reproduction reports from isolated completed runs, never paper files."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd

from .assets import sha
from .reproduction_primary import write_json
from .reproduction_replay import patients, old_anchor_path, PRIMARY


def comparison(rows, study, patient, method, mask, metric, paper, rerun, aggregation, status, **extra):
    paper, rerun = float(paper), float(rerun)
    rows.append({'study':study,'patient':patient,'method':method,'mask':mask,'metric':metric,
                 'paper_value':paper,'rerun_value':rerun,'delta':rerun-paper,'absolute_delta':abs(rerun-paper),
                 'relative_delta':(rerun-paper)/abs(paper) if paper else None,
                 'aggregation':aggregation,'status':status,**extra})


def generate(root):
    out=root/'reproduction';runs=out/'runs';reports=out/'verification';reports.mkdir(exist_ok=True)
    stages={'tier_a':'tier-a-v1','tier_b':'tier-b-v2','tier_c':'tier-c-v3','matched':'matched-v2','training':'primary-20261101-v1'}
    summaries={}
    for name,folder in stages.items():
        directory=runs/folder
        for p in directory.glob('*.json'):
            shutil.copyfile(p,reports/(name+'_'+p.name))
        if (directory/'SUMMARY.json').exists():summaries[name]=json.loads((directory/'SUMMARY.json').read_text())
    for folder in ['tier-b-v2','matched-v2']:
        for p in (runs/folder).glob('*.csv'):shutil.copyfile(p,reports/(folder+'_'+p.name))
    for folder in ['tier-c-v2','matched-v1']:
        for name in ['SUMMARY.json','REGRESSION.json','ARTIFACT_CHECKS.json']:
            p=runs/folder/name
            if p.exists():shutil.copyfile(p,reports/('runtime_sensitivity_'+folder+'_'+name))
    training=json.loads((runs/'primary-20261101-v1/TRAINING_MANIFEST.json').read_text())
    rows=[]
    # Full twenty-estimator20% patient comparison, with all mask outcomes in the result asset.
    for cohort,patient in patients(root):
        actual=pd.read_csv(runs/'tier-a-v1'/patient/'PANEL_METRICS.csv.gz')
        expected=pd.read_csv(old_anchor_path(root,cohort,patient)/'PANEL_METRICS.csv.gz')
        for label,frame in [('a',actual),('b',expected)]:
            subset=frame[frame.experiment.eq('A')&frame.budget.eq(.2)]
            if label=='a':a=subset.groupby('method')[PRIMARY].mean()
            else:b=subset.groupby('method')[PRIMARY].mean()
        for method in a.index:
            for metric in PRIMARY:
                comparison(rows,'tier_a_twenty_estimators',patient,method,'original V1 seeds0-9;20%',metric,b.loc[method,metric],a.loc[method,metric],
                           'equal mean over ten masks within patient','EXACT_REPLAY',reference_kind='frozen journal numerical export')
    b=pd.read_csv(runs/'tier-b-v2/PATIENT_COMPARISON.csv')
    for _,r in b.iterrows():
        for metric in PRIMARY:
            comparison(rows,'tier_b_primary',r.patient,r.method,'original V1 seeds0-9;20%',metric,r[metric+'_paper'],r[metric+'_rerun'],
                       'equal mean over ten masks within patient','APPROXIMATE_RETRAIN',training_seed=r.training_seed,reference_kind='original C9/C10 numerical export')
    panel=pd.read_csv(runs/'tier-b-v2/PANEL.csv.gz')
    original=[]
    for cohort,cp in [('wu','C9'),('bassiouni','C10')]:
        d=pd.read_csv(root/f'checkpoints/legacy/{cp}/{cp}B/{cp}B_PANEL_METRICS.csv')
        d=d[d.regime.eq('V1')&d.strategy.eq('random')&d.budget.eq(.2)&d.control.fillna('').eq('')&d.method.isin(['M0','M3'])].copy()
        d['patient']=d['sample_id'] if 'sample_id' in d else d['patient_id'];d=d[d.patient.isin(panel.patient)]
        original.append(d.rename(columns={'seed':'mask_seed'}))
    reference=pd.concat(original)
    merged=panel.merge(reference[['patient','method','mask_seed']+PRIMARY],on=['patient','method','mask_seed'],suffixes=('_rerun','_paper'),validate='many_to_one')
    for _,r in merged.iterrows():
        for metric in PRIMARY:
            comparison(rows,'tier_b_primary_mask',r.patient,r.method,r['mask'],metric,r[metric+'_paper'],r[metric+'_rerun'],
                       'per-mask panel metric','APPROXIMATE_RETRAIN',training_seed=r.training_seed,reference_kind='original C9/C10 numerical export')
    cohort=pd.read_csv(runs/'tier-b-v2/COHORT_COMPARISON.csv')
    for _,r in cohort.iterrows():
        comparison(rows,'tier_b_primary_cohort',r.cohort,r.method,'V1;20%',r.metric,r.paper_value,r.rerun_value,r.aggregation,r.status,
                   training_seed=r.training_seed,reference_kind='full precision original C9/C10 cohort aggregation',criterion_pass=r.criterion_pass)
    patient=pd.read_csv(runs/'tier-a-v1/PATIENT_METRICS.csv')
    for cohort,paper in [('Wu',.359),('Bassiouni',.416)]:
        d=patient[patient.cohort.eq(cohort)&patient.experiment.eq('A')&patient.budget.eq(.2)&patient.method.eq('idw_k16_p2')]
        comparison(rows,'tier_a_paper_anchor',cohort,'idw_k16_p2','V1;20%','median_pearson',paper,d.median_pearson.median(),
                   'mean masks within patient; median patients','EXACT_REPLAY',reference_kind='typeset three-decimal target; strict full-precision audit separate')
    effects=pd.read_csv(runs/'tier-a-v1/PAIRED_PATIENT_EFFECTS.csv')
    d=effects[effects.cohort.eq('Bassiouni')&effects.budget.eq(.2)&effects.contrast.eq('idw_k16_p2_region_exclusion_vs_distributed')]
    for metric,paper in [('delta_pearson',-.1541),('relative_mse',.45)]:
        comparison(rows,'tier_a_paper_placement','Bassiouni','region_exclusion_vs_distributed','B2 identical H;20%',metric,paper,d[metric].median(),
                   'paired per-patient effect; median22 patients','EXACT_REPLAY',reference_kind='rounded typeset target; strict full-precision audit separate')
    matched=pd.read_csv(runs/'matched-v2/COHORT_EFFECTS.csv').set_index('backbone')
    for backbone,targets in [('resnet18',(.0481,-.1784)),('phikon',(.0075,-.0170))]:
        for metric,paper in zip(['delta_pearson','relative_mse'],targets):
            comparison(rows,'matched_recovered_weights','Bassiouni',backbone+'_GCN_vs_local','historical neural seeds0-2;20% residual',metric,
                       paper,matched.loc[backbone,metric],'equal masks, equal training seeds, paired patient effect; median22 patients','EXACT_REPLAY',
                       reference_kind='rounded typeset target;792 exact full-precision reference rows audited separately')
    pd.DataFrame(rows).to_csv(out/'RESULTS_COMPARISON.csv',index=False)
    # Append replacement identities without changing or erasing the original missing rows.
    inventory=list(csv.DictReader(open(out/'INVENTORY.csv')))
    inventory=[r for r in inventory if r['family']!='retrained_primary']
    cp=root/'checkpoints/reproduced/primary_resnet18/primary-20261101-v1'
    for r in training['heads']:
        inventory.append({'identity':'replacement_primary_seed'+str(r['seed']),'family':'retrained_primary','status':'RETRAINED_REPLACEMENT',
                          'expected_sha256':'','actual_sha256':r['sha256'],'path':str((cp/r['file']).relative_to(root)),
                          'needed_for_current_replay':True,'note':'RETRAINED_NOT_ORIGINAL; fresh identity; source/validation only; original remains missing.'})
    with (out/'INVENTORY.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(inventory[0]),lineterminator='\n');w.writeheader();w.writerows(inventory)
    # Model-specific index, preserving candidate/selected distinctions and shared state references.
    items={}; selected={r['run_id'] for r in json.loads((root/'experiments/journal_revision/backbone_graph_comparison/SELECTED_FIT_MANIFEST.json').read_text())}
    sources=[('minibatch_30',root/'experiments/journal_revision/results/hest_learned/FIT_MANIFEST.json'),
             ('matched_candidates_72',root/'experiments/journal_revision/backbone_graph_comparison/ALL_CANDIDATE_FIT_MANIFEST.json')]
    for study,path in sources:
        for fit in json.loads(path.read_text()):
            checkpoint=fit.get('checkpoint',fit.get('checkpoint_relative_path'));backbone=fit.get('backbone','resnet18')
            d=items.setdefault(checkpoint,{'checkpoint':checkpoint,'sha256':fit.get('checkpoint_sha256'),'backbone':backbone,'origin':'RECOVERED_POST_V1_EXPLORATORY',
                'seed':fit['seed'],'model':fit['model'],'control':fit['control'],'architecture':{'class':'core.mlp_model','input_dim':(768 if backbone=='phikon' else 512)*(2 if fit['model']=='SAGE' else 1),'output_dim':50,**fit['config']},
                'feature_scaler':fit.get('feature_scaler','checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz'),
                'target_scaling':'Direct C0 normalized targets; no additional head-target scaler',
                'selected_epoch':fit.get('selected_epoch',fit.get('best_epoch')),'training':fit.get('training'),
                'fit_records':[],'references':[]})
            d['fit_records'].append({'source':str(path.relative_to(root)),'sha256':sha(path),'record':fit})
            d['references'].append({'study':study,'run_id':fit['run_id'],'matched_selected':study=='matched_candidates_72' and fit['run_id'] in selected})
    for path in sorted((root/'experiments/journal_revision/transformer_egnn/fits').glob('*.json')):
        fit=json.loads(path.read_text());checkpoint=fit['best_checkpoint_relative']
        items[checkpoint]={'checkpoint':checkpoint,'sha256':fit['best_checkpoint_sha256'],'backbone':'resnet18','origin':'RECOVERED_POST_V1_EXPLORATORY',
            'seed':fit['seed'],'model':fit['variant'],'control':'see immutable variant definition','architecture':{'class':'extension_models.ExtensionModel','variant':fit['variant']},
            'feature_scaler':'checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz','target_scaling':'Direct C0 normalized targets; no additional head-target scaler',
            'selected_epoch':fit['selected_epoch'],'training':fit.get('config'),
            'fit_records':[{'source':str(path.relative_to(root)),'sha256':sha(path),'record':fit}],
            'references':[{'study':'advanced_48','run_id':fit['variant']+'_seed'+str(fit['seed']),'matched_selected':False}]}
    checks=json.loads((runs/'tier-c-v3/ARTIFACT_CHECKS.json').read_text())
    for path,d in items.items():
        d['bytes']=(root/path).stat().st_size
        if sha(root/path)!=d['sha256']:raise ValueError('Recovered model changed after verification')
        d['verification']=[r for r in checks if r['checkpoint']==path]
        d['distribution']='DRAFT_RIGHTS_REVIEW' if d['backbone']=='phikon' else 'CC-BY-NC-SA-4.0'
    index={'schema_version':1,'genes':training['genes'],'fit_samples':['NCBI783','NCBI785'],'validation_sample':'TENX95',
           'models':list(items.values()),'logical_references':150,'unique_best_states':len(items),
           'pretrained_encoders_included':False,'original_primary_weights_recovered':False}
    write_json(out/'MODEL_INDEX.json',index)
    write_json(out/'COMPUTE_SUMMARY.json',{'training_elapsed_seconds':training['elapsed_seconds'],'training_peak_rss_bytes':training['peak_rss_bytes'],
               'environment':training['environment'],'stages':summaries,'new_training_fits':3,'public_research_dataset_bytes_downloaded_this_task':0,
               'paid_computing_services_used':False,'note':'Research inputs were recovered caches; fresh license files and the clean-checkout encoder download are recorded separately.'})
    return {'comparison_rows':len(rows),'models':len(items),'public_recovered_heads':sum(d['backbone']=='resnet18' for d in items.values()),
            'draft_phikon_heads':sum(d['backbone']=='phikon' for d in items.values()),'stages':summaries}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);a=p.parse_args()
    print(json.dumps(generate(a.root.resolve()),indent=2))
