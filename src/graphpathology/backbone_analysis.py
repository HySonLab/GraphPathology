"""Patient-paired graph/local/anchor contrasts with explicit sensitivity analysis."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .assets import sha
from .backbone_data import STUDY,config,record,write
from .extension_evaluate import cohort
from .journal_methods import patient_bootstrap

METRICS=['median_pearson','mean_mse','mean_mae']
DIMS=['backbone','cohort','patient','model','control','route','design','condition','budget']
PAIR=['cohort','patient','route','design','condition','budget']
RUN=['patient','design','condition','budget','mask_seed','geometry','training_seed','mask_sha256']


def patient_average(frame, values=METRICS):
    """Equal masks within geometry, then geometries, then training seeds."""
    a=frame.groupby(DIMS+['training_seed','geometry'],dropna=False)[values].mean().reset_index()
    b=a.groupby(DIMS+['training_seed'],dropna=False)[values].mean().reset_index()
    return b.groupby(DIMS,dropna=False)[values].mean().reset_index()


def paired_effects(candidate,reference):
    if candidate.duplicated(PAIR).any() or reference.duplicated(PAIR).any():raise ValueError('Nonunique patient pairing')
    joined=candidate.merge(reference,on=PAIR,suffixes=('_candidate','_reference'),validate='one_to_one',how='outer',indicator=True)
    if not (joined['_merge']=='both').all():raise ValueError('Unmatched patient/support conditions')
    joined['delta_pearson']=joined.median_pearson_candidate-joined.median_pearson_reference
    joined['delta_mae']=joined.mean_mae_candidate-joined.mean_mae_reference
    joined['relative_mse']=joined.mean_mse_candidate/joined.mean_mse_reference-1
    return joined


def common_valid_runs(candidate,reference,baseline=False):
    keys=[k for k in RUN if not(baseline and k=='training_seed')]+['gene']
    if candidate.duplicated(RUN+['gene']).any() or reference.duplicated(keys).any():raise ValueError('Nonunique gene/support pairing')
    merged=candidate.merge(reference[keys+['pearson']],on=keys,suffixes=('_candidate','_reference'),validate='many_to_one' if baseline else 'one_to_one',how='left',indicator=True)
    if not (merged['_merge']=='both').all():raise ValueError('Gene/support mismatch in common-valid sensitivity')
    merged['valid']=np.isfinite(merged.pearson_candidate)&np.isfinite(merged.pearson_reference)
    rows=[]
    for identity,g in merged.groupby(RUN,dropna=False):
        good=g[g.valid]
        rows.append(dict(zip(RUN,identity))|{'common_valid_genes':len(good),'candidate_defined_genes':int(np.isfinite(g.pearson_candidate).sum()),'reference_defined_genes':int(np.isfinite(g.pearson_reference).sum()),'common_delta_pearson':float(good.pearson_candidate.median()-good.pearson_reference.median()) if len(good) else np.nan})
    return pd.DataFrame(rows)


def manifests(root):
    base=root/f'experiments/journal_revision/{STUDY}'
    reused=json.loads((base/'REUSED_EVALUATION_MANIFEST.json').read_text())
    new_path=base/'EVALUATION_MANIFEST.json'
    new=json.loads(new_path.read_text()) if new_path.exists() else []
    return reused+new


def bundle_paths(root,item):
    if item['status']=='REUSED_VERIFIED':return root/item['panel'],root/item['genes']
    p=root/f'experiments/journal_revision/{STUDY}/evaluation'/item['patient']/item['run_id']
    for name,digest in item['output_hashes'].items():assert sha(p/name)==digest
    return p/'PANEL.csv.gz',p/'GENES.csv.gz'


def collect(root,items):
    frames=[]
    for item in items:
        panel,_=bundle_paths(root,item);f=pd.read_csv(panel)
        for k in ['backbone','model','control']:f[k]=item[k]
        f['training_seed']=item['training_seed'];f['undefined_pearson_genes']=f.n_genes-f.n_defined_pearson;frames.append(f)
    for patient in sorted({r['patient'] for r in items}):
        p=root/'experiments/journal_revision/transformer_egnn/evaluation'/patient/'MOLECULAR.csv'
        f=pd.read_csv(p);f['backbone']='molecular';f['model']=f.variant;f['control']='none';f['training_seed']=-1;f['cohort']=cohort(patient)
        f['undefined_pearson_genes']=f.n_genes-f.n_defined_pearson;frames.append(f)
    return pd.concat(frames,ignore_index=True)


def references(frame,model,control,route):
    options=[('local','local',route)] if model!='local' else []
    if model in ['GCN','SAGE'] and control=='spatial':options += [('MEAN','spatial',route),(model,'shuffled',route)]
    if route!='image':options += [('idw_k16_p2','none','anchor_only'),('uniform_k16','none','anchor_only')]
    return options


def contrasts(patient):
    rows=[]
    for (backbone,model,control,route),candidate in patient[patient.backbone!='molecular'].groupby(['backbone','model','control','route']):
        for reference_model,reference_control,reference_route in references(patient,model,control,route):
            ref=patient[(patient.backbone==('molecular' if reference_control=='none' else backbone))&(patient.model==reference_model)&(patient.control==reference_control)&(patient.route==reference_route)].copy()
            ref['route']=route
            paired=paired_effects(candidate,ref)
            paired['backbone']=backbone;paired['candidate']=model+'_'+control;paired['reference']=reference_model+'_'+reference_control
            rows.append(paired[['backbone']+PAIR+['candidate','reference','delta_pearson','relative_mse','delta_mae']])
    return pd.concat(rows,ignore_index=True)


def summaries(effects,cfg):
    rows=[];keys=['backbone','cohort','design','condition','budget','route','candidate','reference']
    for identity,g in effects.groupby(keys,dropna=False):
        for metric in ['delta_pearson','relative_mse','delta_mae']:
            values=g[metric].to_numpy(float);finite=values[np.isfinite(values)]
            if identity[1]=='Bassiouni':stats=patient_bootstrap(values,'positive' if metric=='delta_pearson' else 'negative',cfg['bootstrap_resamples'],cfg['bootstrap_seed'])
            else:stats={'n_patients':len(finite),'median_effect':float(np.median(finite)) if len(finite) else np.nan,'ci_low':np.nan,'ci_high':np.nan,'bootstrap_resamples':0,'bootstrap_seed':cfg['bootstrap_seed']}
            stats['n_improved']=int(np.sum(finite>0 if metric=='delta_pearson' else finite<0))
            rows.append(dict(zip(keys,identity))|{'effect':metric,**stats})
    return pd.DataFrame(rows)


def sensitivity(root,items):
    chunks=[];meta=[]
    for patient in sorted({r['patient'] for r in items}):
        molecular=pd.read_csv(root/'experiments/journal_revision/transformer_egnn/evaluation'/patient/'MOLECULAR_GENES.csv.gz')
        for backbone in sorted({r['backbone'] for r in items if r['patient']==patient}):
            selected=[r for r in items if r['patient']==patient and r['backbone']==backbone];frames=[]
            for item in selected:
                _,gene_path=bundle_paths(root,item);g=pd.read_csv(gene_path);g['model']=item['model'];g['control']=item['control'];g['training_seed']=item['training_seed'];frames.append(g)
            genes=pd.concat(frames,ignore_index=True)
            for (model,control,route),candidate in genes.groupby(['model','control','route']):
                for rm,rc,rr in references(genes,model,control,route):
                    baseline=rc=='none'
                    ref=molecular[molecular.model_id==rm] if baseline else genes[(genes.model==rm)&(genes.control==rc)&(genes.route==rr)]
                    result=common_valid_runs(candidate,ref,baseline)
                    result['backbone']=backbone;result['cohort']=cohort(patient);result['candidate']=model+'_'+control;result['reference']=rm+'_'+rc;result['route']=route
                    chunks.append(result)
            meta.append({'patient':patient,'backbone':backbone,'configurations':len(selected)})
        print('COMMON VALID SENSITIVITY',patient,flush=True)
    return pd.concat(chunks,ignore_index=True),meta


def analyze(root,include_sensitivity=True):
    cfg=config(root);items=manifests(root);out=root/f'experiments/journal_revision/{STUDY}/summary';out.mkdir(parents=True,exist_ok=True)
    frames=collect(root,items);patients=patient_average(frames);effects=contrasts(patients);summary=summaries(effects,cfg)
    frames.to_csv(out/'RUN_METRICS.csv.gz',index=False);patients.to_csv(out/'PATIENT_METRICS.csv',index=False);effects.to_csv(out/'PAIRED_PATIENT_EFFECTS.csv',index=False);summary.to_csv(out/'PAIRED_SUMMARY.csv',index=False)
    undefined=frames.groupby(['backbone','cohort','model','control','route'],dropna=False).agg(run_rows=('n_defined_pearson','size'),minimum_defined_genes=('n_defined_pearson','min'),maximum_defined_genes=('n_defined_pearson','max'),undefined_gene_instances=('undefined_pearson_genes','sum')).reset_index();undefined.to_csv(out/'DEFINED_GENE_COUNTS.csv',index=False)
    if include_sensitivity:
        runs,meta=sensitivity(root,items);runs.to_csv(out/'COMMON_VALID_RUN_EFFECTS.csv.gz',index=False)
        dims=['backbone','cohort','patient','candidate','reference','route','design','condition','budget']
        a=runs.groupby(dims+['training_seed','geometry'],dropna=False)[['common_delta_pearson','common_valid_genes']].mean().reset_index()
        b=a.groupby(dims+['training_seed'],dropna=False)[['common_delta_pearson','common_valid_genes']].mean().reset_index()
        b.groupby(dims,dropna=False)[['common_delta_pearson','common_valid_genes']].mean().reset_index().to_csv(out/'COMMON_VALID_PATIENT_EFFECTS.csv',index=False)
    report=record(root,status='COMPLETE' if sum(r['backbone']=='phikon' for r in items)==468 else 'REUSED_VERIFIED_RESNET_ONLY_PHIKON_PENDING',reused_bundles=sum(r['status']=='REUSED_VERIFIED' for r in items),new_bundles=sum(r['status']=='COMPLETE' for r in items),patient_effect_rows=len(effects),source_pointers_only_for_reused_gene_arrays=True,output_hashes={f.name:sha(f) for f in out.iterdir() if f.is_file()})
    write(root/f'experiments/journal_revision/{STUDY}/ANALYSIS_STATUS.json',report)
    return report


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--without-sensitivity',action='store_true');a=p.parse_args();print(json.dumps(analyze(a.root.resolve(),not a.without_sensitivity),indent=2))
