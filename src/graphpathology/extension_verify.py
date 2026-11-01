"""Offline integrity/status verification for the separate extension evidence."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .assets import sha
from .extension_train import metadata

def run(root,require_complete=True):
 base=root/'experiments/journal_revision/transformer_egnn';snapshot=json.loads((base/'BASELINE.json').read_text())
 # Accept the frozen snapshot schema explicitly, never infer a new expected hash from current files.
 original=snapshot['file_hashes']
 if isinstance(original,list):original={r['path']:r['sha256'] for r in original}
 if not original:raise ValueError('Missing original preservation manifest')
 differences=[str(p) for p,digest in original.items() if not (root/p).exists() or sha(root/p)!=digest]
 if differences:raise ValueError('Protected evidence changed: '+str(differences))
 failures=[];checked=0
 for p in base.rglob('COMPLETE.json'):
  r=json.loads(p.read_text())
  if 'output_hashes' not in r:raise ValueError('Unhashed completion record '+str(p))
  for name,digest in r['output_hashes'].items():
   if sha(p.parent/name)!=digest:failures.append(str(p.parent/name))
  checked+=1
 if failures:raise ValueError('Modified outputs: '+str(failures))
 diag=json.loads((base/'diagnostics/STATUS.json').read_text());assert diag['status']=='COMPLETE'
 if require_complete:
  summary=json.loads((base/'summary/STATUS.json').read_text());assert summary['status']=='COMPLETE' and summary['models']==78 and summary['new_fits']==48 and summary['patients']==26
  assert len(list((base/'bassiouni_anchors').glob('*/COMPLETE.json')))==22
 report=metadata(root,status='COMPLETE' if require_complete else 'IN_PROGRESS',protected_files_unchanged=len(original),completed_bundles_hash_verified=checked,geometric_checks=diag['geometric_checks'],raw_datasets_and_weights_in_git=False)
 (base/'VALIDATION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n');return report
