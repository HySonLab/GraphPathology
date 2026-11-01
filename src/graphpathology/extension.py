"""Explicit post_v1_exploratory extension commands. Imports/help/tests never acquire or fit."""
import argparse,json,traceback
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);s=p.add_subparsers(dest='action',required=True)
 s.add_parser('preflight');s.add_parser('prepare-wu');s.add_parser('prepare-bassiouni');s.add_parser('acquire-bassiouni');q=s.add_parser('train');q.add_argument('--variant',required=True);q.add_argument('--seed',type=int,choices=[0,1,2],required=True);q.add_argument('--reproduction-output',type=Path)
 q=s.add_parser('train-all');q.add_argument('--reproduction-output',type=Path);q=s.add_parser('evaluate');q.add_argument('--patient',action='append');q.add_argument('--threads',type=int,default=4);s.add_parser('diagnostics');s.add_parser('benchmark');q=s.add_parser('paper');q.add_argument('--build-output',type=Path);s.add_parser('report');q=s.add_parser('verify');q.add_argument('--committed-only',action='store_true');q=s.add_parser('bundle');q.add_argument('--output',type=Path,required=True);q=s.add_parser('source-package');q.add_argument('--output',type=Path,required=True);q.add_argument('--check-directory',type=Path,required=True)
 a=p.parse_args();r=a.root.resolve()
 if a.action=='preflight':
  from .extension_train import preflight
  print(json.dumps(preflight(r),indent=2))
 elif a.action=='acquire-bassiouni':
  from .extension_geo import acquire
  acquire(r)
 elif a.action=='prepare-bassiouni':
  from .extension_data import prepare_bassiouni
  prepare_bassiouni(r)
 elif a.action=='prepare-wu':
  from .extension_data import prepare_wu
  prepare_wu(r)
 elif a.action in ['train','train-all']:
  from .extension_models import VARIANTS
  from .extension_train import train,metadata
  if a.reproduction_output:
   from .extension_reproduction import staging_root
   r=staging_root(r,a.reproduction_output)
  variants=[a.variant] if a.action=='train' else list(VARIANTS);seeds=[a.seed] if a.action=='train' else [0,1,2]
  for seed in seeds:
   for variant in variants:
    if variant not in VARIANTS:raise ValueError('Unknown frozen variant')
    try:train(r,variant,seed)
    except Exception as exc:
     path=r/'experiments/journal_revision/transformer_egnn/fits'/f'{variant}_seed{seed}_FAILED.json';path.parent.mkdir(exist_ok=True)
     if not path.exists():path.write_text(json.dumps(metadata(r,status='FAILED_VALIDATION',variant=variant,seed=seed,error=str(exc),traceback=traceback.format_exc()),indent=2)+'\n')
     if a.action=='train':raise
 elif a.action=='source-package':
  from .extension_source import run
  print(json.dumps(run(r,a.output,a.check_directory),indent=2))
 elif a.action=='bundle':
  from .extension_bundle import run
  result=run(r,a.output);print(json.dumps({'status':result['status'],'bundles':result['bundles']},indent=2))
 elif a.action=='paper':
  from .extension_paper import build
  print(json.dumps(build(r,output=a.build_output),indent=2))
 elif a.action in ['evaluate','diagnostics','benchmark','report','verify']:
  from importlib import import_module
  module=import_module('.extension_'+{'evaluate':'evaluate','diagnostics':'diagnostics','benchmark':'benchmark','report':'report','verify':'verify'}[a.action],package='graphpathology')
  print(json.dumps(module.run(r,patients=a.patient,threads=a.threads) if a.action=='evaluate' else module.run(r,committed_only=a.committed_only) if a.action=='verify' else module.run(r),indent=2))

if __name__=='__main__':main()
