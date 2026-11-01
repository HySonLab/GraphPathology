"""Verify the selected two-rate study; private arrays are explicitly skipped in CI."""
import json
from pathlib import Path
import pandas as pd
from .assets import sha
from .backbone_data import STUDY,config,record,write
from .backbone_verify import verify as verify_initial
from .backbone_refinement import validation_selection


def verify(root,committed_only=False,require_complete=True):
    base=root/f'experiments/journal_revision/{STUDY}';cfg=config(root)
    prior=verify_initial(root,committed_only,require_complete=False)
    candidates=json.loads((base/'ALL_CANDIDATE_FIT_MANIFEST.json').read_text());selected=json.loads((base/'SELECTED_FIT_MANIFEST.json').read_text());assert len(candidates)==72 and len(selected)==36
    expected,selection=validation_selection(candidates)
    signature=lambda values:sorted((r['backbone'],r['model'],r['control'],r['seed'],r['learning_rate'],r['checkpoint_sha256']) for r in values)
    assert signature(expected)==signature(selected)
    assert selection==json.loads((base/'VALIDATION_SELECTION.json').read_text())['records']
    states=skipped_states=0
    for r in candidates:
        checkpoint=root/r['checkpoint']
        if checkpoint.exists():assert sha(checkpoint)==r['checkpoint_sha256'];states+=1
        elif committed_only:skipped_states+=1
        else:raise FileNotFoundError(checkpoint)
        assert r['training']['max_epochs']==80 and r['training']['batch_size']==512
        assert r['learning_rate'] in [.001,.0001]
    path=base/'PRIMARY_EVALUATION_MANIFEST.json';evaluations=json.loads(path.read_text()) if path.exists() else []
    if require_complete:assert len(evaluations)==936
    lookup={r['run_id']:r for r in selected};genes_checked=gene_skips=0
    for r in evaluations:
        assert r['run_id'] in lookup
        item=lookup[r['run_id']]
        assert r['backbone']==item['backbone'] and r['training_seed']==item['seed'] and r['learning_rate']==item['learning_rate']
        assert r['selected_fit_manifest_sha256']==sha(base/'SELECTED_FIT_MANIFEST.json')
        for name in ['panel','genes']:
            p=root/r[name]
            if not p.exists() and committed_only and name=='genes':gene_skips+=1;continue
            assert sha(p)==r[name+'_sha256']
            if name=='genes':genes_checked+=1
        if r['status']=='COMPLETE':
            for name,digest in r['output_hashes'].items():
                p=(root/r['panel']).parent/name
                if not p.exists() and committed_only and name=='GENES.csv.gz':continue
                assert sha(p)==digest
            assert r['calibration_applications']==(0 if r['patient']=='TENX99' else 1)
            assert r['exact_label_corruption_checks']==4*r['mask_sets'] and r['image_anchor_independence_checks']>0
            assert r['checkpoint_sha256']==item['checkpoint_sha256']
        else:assert r['status']=='REUSED_VERIFIED' and r['prior_exact_label_corruption_checks']>0
    expected_configs={(c['model'],c['control'],s) for c in cfg['configurations'] for s in cfg['seeds']}
    if require_complete:
        for patient in ['TENX99']+cfg['external']:
            for backbone in cfg['backbones']:assert {(r['model'],r['control'],r['training_seed']) for r in evaluations if r['patient']==patient and r['backbone']==backbone}==expected_configs
        analysis=json.loads((base/'ANALYSIS_STATUS.json').read_text());assert analysis['status']=='COMPLETE'
        for name,digest in analysis['output_hashes'].items():assert sha(base/'summary'/name)==digest
        effect=pd.read_csv(base/'summary/PAIRED_PATIENT_EFFECTS.csv');keys=['backbone','cohort','patient','route','design','condition','budget','candidate','reference'];assert not effect.duplicated(keys).any()
        counts=effect.groupby([k for k in keys if k!='patient']).size().reset_index(name='n');assert (counts.n==counts.cohort.map({'Development':1,'Wu':3,'Bassiouni':22})).all()
        sensitivity=pd.read_csv(base/'summary/COMMON_VALID_RUN_EFFECTS.csv.gz');assert sensitivity.common_valid_genes.between(0,50).all()
    return record(root,status='COMPLETE' if require_complete else 'IN_PROGRESS_PRESERVATION_PASS',initial_evidence_preservation=prior,candidate_fits=72,validation_selected_fits=36,checkpoint_states_sha256_verified=states,private_checkpoint_states_absent=skipped_states,selected_evaluation_bundles=len(evaluations),new_selected_bundles=sum(r['status']=='COMPLETE' for r in evaluations),reused_selected_bundles=sum(r['status']=='REUSED_VERIFIED' for r in evaluations),gene_export_arrays_sha256_verified=genes_checked,private_gene_arrays_absent=gene_skips,selection_recomputed_from_validation_only=True,all_three_seeds_per_selected_configuration=True)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--committed-only',action='store_true');p.add_argument('--allow-in-progress',action='store_true');p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify(a.root.resolve(),a.committed_only,not a.allow_in_progress);write(a.output,r);print(json.dumps(r,indent=2))
