"""Offline integrity/status verification for the separate extension evidence."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .assets import sha
from .extension_train import metadata

def run(root,require_complete=True,committed_only=False):
 base=root/'experiments/journal_revision/transformer_egnn';snapshot=json.loads((base/'BASELINE.json').read_text())
 # Accept the frozen snapshot schema explicitly, never infer a new expected hash from current files.
 original=snapshot['file_hashes']
 if isinstance(original,list):original={r['path']:r['sha256'] for r in original}
 if not original:raise ValueError('Missing original preservation manifest')
 skipped_checkpoints=[];differences=[]
 for p,digest in original.items():
  if committed_only and not (root/p).exists() and p.startswith('checkpoints/journal_revision/hest_learned/'):
   skipped_checkpoints.append(p);continue
  if not (root/p).exists() or sha(root/p)!=digest:differences.append(p)
 if differences:raise ValueError('Protected evidence changed: '+str(differences))
 failures=[];checked=0;skipped_genes=[]
 release=json.loads((base/'RELEASE_READY_MANIFEST.json').read_text()) if (base/'RELEASE_READY_MANIFEST.json').exists() else {}
 bundled={r['path']:r['sha256'] for r in release.get('full_gene_files',[])}
 for p in base.rglob('COMPLETE.json'):
  r=json.loads(p.read_text())
  if 'output_hashes' not in r:raise ValueError('Unhashed completion record '+str(p))
  for name,digest in r['output_hashes'].items():
   path=p.parent/name
   if committed_only and not path.exists() and str(path.relative_to(root)) in bundled:
    assert bundled[str(path.relative_to(root))]==digest;skipped_genes.append(str(path.relative_to(root)));continue
   if sha(path)!=digest:failures.append(str(path))
  checked+=1
 if failures:raise ValueError('Modified outputs: '+str(failures))
 diag=json.loads((base/'diagnostics/STATUS.json').read_text());assert diag['status']=='COMPLETE'
 for name,digest in diag['output_hashes'].items():assert sha(base/'diagnostics'/name)==digest
 for p in (base/'evaluation').glob('*/MOLECULAR_CHECKS.json'):
  r=json.loads(p.read_text());assert r['status']=='COMPLETE'
  for name,digest in r['output_hashes'].items():
   path=p.parent/name
   if committed_only and not path.exists() and str(path.relative_to(root)) in bundled:
    assert bundled[str(path.relative_to(root))]==digest;skipped_genes.append(str(path.relative_to(root)));continue
   assert sha(path)==digest
 for p in (base/'fits').glob('*.json'):
  r=json.loads(p.read_text());assert r['status']=='COMPLETE' and len(r['history'])==r['epochs_ran']<=80 and r['optimizer_state_included']
  assert r['validation_mse']==r['history'][r['selected_epoch']-1]['validation_mse']
  assert sha(root/'src/graphpathology/extension_models.py')==r['code_hashes']['src/graphpathology/extension_models.py']
  if not committed_only:assert sha(root/r['best_checkpoint_relative'])==r['best_checkpoint_sha256']
 if require_complete:
  summary=json.loads((base/'summary/STATUS.json').read_text());assert summary['status']=='COMPLETE' and summary['models']==78 and summary['new_fits']==48 and summary['patients']==26
  assert len(list((base/'bassiouni_anchors').glob('*/COMPLETE.json')))==22
 report=metadata(root,status=('COMPLETE_COMMITTED_EXPORTS' if committed_only else 'COMPLETE') if require_complete else 'IN_PROGRESS',skipped_author_local_old_checkpoints=skipped_checkpoints,full_gene_bundle_files_not_installed=skipped_genes,protected_files_unchanged=len(original),completed_bundles_hash_verified=checked,geometric_checks=diag['geometric_checks'],raw_datasets_and_weights_in_git=False)
 target=base/'VALIDATION_REPORT.json'
 if committed_only:
  (root/'outputs').mkdir(exist_ok=True);target=root/'outputs/extension-validation-committed.json'
 target.write_text(json.dumps(report,indent=2)+'\n');return report


def audit_bassiouni(root):
 """All historical random molecular references, not only designated seed0 replays."""
 base=root/'experiments/journal_revision/transformer_egnn';source=base/'bassiouni_anchors';old=root/'checkpoints/legacy/C10/C10B';errors={};counts={}
 for scope,newname,oldname,metrics,expected_rows in [
  ('panel','PANEL_METRICS.csv.gz','C10B_PANEL_METRICS.csv',['mean_pearson','median_pearson','mean_mse','median_mse','mean_mae','median_mae'],2640),
  ('gene','GENE_METRICS.csv.gz','C10B_GENE_METRICS.csv.gz',['pearson','mse','mae'],132000)]:
  actual=pd.concat([pd.read_csv(p) for p in sorted(source.glob('*/'+newname))],ignore_index=True);actual=actual[((actual.regime=='V1')&(actual.method=='idw_k16_p2'))|((actual.regime=='V2')&(actual.method=='uniform_k16'))]
  expected=pd.read_csv(old/oldname);expected=expected[(expected.method=='M1')&(expected.strategy=='random')&expected.control.fillna('').eq('')].rename(columns={'patient_id':'patient'});keys=['patient','regime','geometry','seed','budget']+(['gene'] if scope=='gene' else [])
  merged=actual[keys+metrics].merge(expected[keys+metrics],on=keys,suffixes=('_actual','_expected'),validate='one_to_one',how='outer',indicator=True);assert len(merged)==expected_rows and merged['_merge'].eq('both').all();counts[scope]=len(merged);errors[scope]={}
  for metric in metrics:
   a=merged[metric+'_actual'];b=merged[metric+'_expected'];np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-10,equal_nan=True);errors[scope][metric]=float(np.nanmax(abs(a-b)))
 result={'status':'PASS_PRIMARY_ALL_ORIGINAL_RANDOM_REFERENCES','run_rows':counts['panel'],'gene_rows':counts['gene'],'primary_gene_values_checked':counts['gene']*3,'max_abs_deltas':errors,'rtol':1e-10,'atol':1e-10,'secondary_spearman':'Retained and flagged near-tie sensitive; no exactness claim','original_panel_sha256':sha(old/'C10B_PANEL_METRICS.csv'),'original_gene_sha256':sha(old/'C10B_GENE_METRICS.csv.gz')}
 (base/'BASSIOUNI_ALL_REFERENCE_REGRESSION.json').write_text(json.dumps(result,indent=2)+'\n');return result
