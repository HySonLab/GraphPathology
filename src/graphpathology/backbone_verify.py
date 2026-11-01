"""Preservation and comparison contracts, with explicit data-free checkout limits."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .assets import sha
from .backbone_data import STUDY,config,record,write


def verify(root,committed_only=False,require_complete=True):
    base=root/f'experiments/journal_revision/{STUDY}';cfg=config(root)
    protected=json.loads((base/'BASELINE.json').read_text())['protected_files']
    released=json.loads((root/'experiments/journal_revision/transformer_egnn/RELEASE_READY_MANIFEST.json').read_text())
    optional={r['path']:r['sha256'] for r in released['full_gene_files']}
    checked=skipped=0
    for name,digest in protected.items():
        path=root/name
        if committed_only and not path.exists() and optional.get(name)==digest:skipped+=1;continue
        assert path.is_file() and sha(path)==digest,name
        checked+=1
    fits=json.loads((base/'REUSED_FIT_MANIFEST.json').read_text());new_path=base/'FIT_MANIFEST.json'
    new=json.loads(new_path.read_text()) if new_path.exists() else []
    if require_complete:assert len(new)==18
    for r in fits+new:
        checkpoint=root/r['checkpoint']
        if checkpoint.exists():assert sha(checkpoint)==r['checkpoint_sha256']
        elif not committed_only:raise FileNotFoundError(checkpoint)
    expected={(c['model'],c['control'],s) for c in cfg['configurations'] for s in cfg['seeds']}
    assert {(r['model'],r['control'],r['seed']) for r in fits}==expected
    if new:
        assert {(r['model'],r['control'],r['seed']) for r in new}==expected
        for r in new:
            assert r['source_fit']==cfg['source_fit'] and r['validation']==cfg['validation']
            assert r['training']==cfg['training'] and r['serialization_exact']
            old=next(f for f in fits if (f['model'],f['control'],f['seed'])==(r['model'],r['control'],r['seed']))
            assert old['training']==r['training']
            for g in r['graphs'].values():
                if 'changed_edge_fraction' in g:assert g['changed_edge_fraction']>.5
    reused=json.loads((base/'REUSED_EVALUATION_MANIFEST.json').read_text());path=base/'EVALUATION_MANIFEST.json'
    evaluations=json.loads(path.read_text()) if path.exists() else []
    assert len(reused)==468
    if require_complete:assert len(evaluations)==468
    gene_arrays_checked=gene_arrays_absent=0
    for r in reused:
        for name in ['panel','genes']:
            p=root/r[name]
            if not p.exists() and committed_only and name=='genes':gene_arrays_absent+=1;continue
            assert sha(p)==r[name+'_sha256']
            if name=='genes':gene_arrays_checked+=1
        assert r['prior_exact_label_corruption_checks']>0
        if r['control']=='shuffled':assert r['graph_check']['changed_edge_fraction']>.5
    for r in evaluations:
        folder=base/'evaluation'/r['patient']/r['run_id']
        for name,digest in r['output_hashes'].items():
            p=folder/name
            if not p.exists() and committed_only and name=='GENES.csv.gz':gene_arrays_absent+=1;continue
            assert sha(p)==digest
            if name=='GENES.csv.gz':gene_arrays_checked+=1
        assert r['calibration_applications']==(0 if r['patient']=='TENX99' else 1)
        assert r['exact_label_corruption_checks']==4*r['mask_sets'] and r['image_anchor_independence_checks']>0
        assert r['checkpoint_sha256']==next(f['checkpoint_sha256'] for f in new if f['run_id']==r['run_id'])
    for patient in ['TENX99']+cfg['external']:
        assert {(r['model'],r['control'],r['training_seed']) for r in reused if r['patient']==patient}==expected
        if require_complete:assert {(r['model'],r['control'],r['training_seed']) for r in evaluations if r['patient']==patient}==expected
        support=json.loads((base/'supports'/(patient+'.json')).read_text())
        assert len(support['masks'])==(27 if patient=='TENX99' else 36)
        for m in support['masks']:assert m['n_observed']>0 and m['n_scored']>0
    for backbone,width in [('resnet18',512),('phikon',768)]:
        h=cfg['hidden'][backbone];local=(width+1)*h['local']+(h['local']+1)*50;sage=(2*width+1)*h['SAGE']+(h['SAGE']+1)*50
        assert abs(local-sage)/local<.001
    if require_complete:
        analysis=json.loads((base/'ANALYSIS_STATUS.json').read_text());assert analysis['status']=='COMPLETE'
        for name,digest in analysis['output_hashes'].items():assert sha(base/'summary'/name)==digest
        patient=pd.read_csv(base/'summary/PATIENT_METRICS.csv');effect=pd.read_csv(base/'summary/PAIRED_PATIENT_EFFECTS.csv')
        assert set(patient.backbone)=={'resnet18','phikon','molecular'}
        assert set(effect.backbone)=={'resnet18','phikon'}
        assert not effect.duplicated(['backbone','cohort','patient','route','design','condition','budget','candidate','reference']).any()
        assert (effect.groupby(['backbone','cohort','route','design','condition','budget','candidate','reference']).size().reset_index(name='n').assign(expected=lambda x:x.cohort.map({'Development':1,'Wu':3,'Bassiouni':22})).eval('n==expected')).all()
    result=record(root,status='COMPLETE' if require_complete else 'REUSED_VERIFIED_WITH_NEW_STUDY_IN_PROGRESS',protected_files_unchanged=checked,author_local_prior_gene_files_skipped=skipped,reused_fits_verified=len(fits),new_fits_verified=len(new),reused_patient_model_bundles=len(reused),new_patient_model_bundles=len(evaluations),gene_export_arrays_checked=gene_arrays_checked,gene_export_arrays_absent_in_checkout=gene_arrays_absent,private_state_verification_skipped=committed_only,all_planned_seeds_retained=True,head_capacity_difference_below_point_one_percent=True)
    return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--committed-only',action='store_true');p.add_argument('--allow-in-progress',action='store_true');p.add_argument('--output',type=Path,required=True);a=p.parse_args();result=verify(a.root.resolve(),a.committed_only,not a.allow_in_progress);write(a.output,result);print(json.dumps(result,indent=2))
