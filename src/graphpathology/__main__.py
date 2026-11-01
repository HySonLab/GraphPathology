"""Offline, explicit entrypoints. No scientific work on import or --help."""
import argparse, json, sys
from pathlib import Path

def main(argv=None):
    parser=argparse.ArgumentParser(description='GraphPathology maintained reproduction tools')
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2],help='checkout root (pass explicitly outside editable installations)')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('doctor',help='Report installed versions and required asset availability')
    p.add_argument('--require-models',action='store_true')
    p=sub.add_parser('verify-archive',help='Read-only byte and historical prefix checks')
    p.add_argument('--historical',action='store_true');p.add_argument('--report',type=Path)
    p=sub.add_parser('resolve',help='Resolve an original path through the migration manifest');p.add_argument('path')
    p=sub.add_parser('smoke',help='Tiny CPU graph and external fixture routes (synthetic heads)');p.add_argument('--output',required=True);p.add_argument('--dry-run',action='store_true')
    p=sub.add_parser('publication',help='Validate current 29 targets or regenerate available exported-data targets')
    p.add_argument('--target',default='all');p.add_argument('--output');p.add_argument('--dry-run',action='store_true');p.add_argument('--available-only',action='store_true')
    p=sub.add_parser('build-documents',help='Build supplied revised LaTeX sources if available')
    p.add_argument('--output');p.add_argument('--dry-run',action='store_true')
    p=sub.add_parser('infer',help='Bounded one-section replay from aligned embeddings and observed targets')
    p.add_argument('--config',required=True,type=Path);p.add_argument('--dry-run',action='store_true')
    args=parser.parse_args(argv);root=args.root.resolve()
    try:
        from . import assets
        if args.command=='resolve':print(assets.resolve(root,args.path));return 0
        if args.command=='verify-archive':
            result=assets.verify(root,args.historical)
            if args.report:
                if args.report.exists():raise FileExistsError('Refusing to overwrite report')
                resolved=args.report.resolve()
                if root in resolved.parents and resolved.relative_to(root).parts[0] not in ['outputs','work']:raise ValueError('Audit reports must use outputs/ or work/')
                resolved.parent.mkdir(parents=True,exist_ok=True);resolved.write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps({k:v for k,v in result.items() if k!='historical_details'},indent=2))
            return int(bool(result['migration']['changed'] or result['migration']['missing'] or result.get('historical_new_issues')))
        if args.command=='doctor':
            from importlib.metadata import version,PackageNotFoundError
            versions={}
            for n in ['graphpathology','numpy','scipy','pandas','torch','torchvision','pytest','numba','matplotlib']:
                try:versions[n]=version(n)
                except PackageNotFoundError:versions[n]='not installed'
            rows=assets.inventory(root);models=[{**x,'available':(root/x['path']).is_file()} for x in rows if x['kind'] in ['learned_state','fitted_pca']]
            primary=[x for x in models if Path(x['path']).name.startswith('C4A_primary_resnet18_reference_primary_S2-50_D1_h256d0.1_s')]
            print(json.dumps({'versions':versions,'recorded_assets':len(rows),'fitted_states':len(models),'installed_fitted_states':sum(x['available'] for x in models),'primary_heads':primary,'publication':str(root/'docs/REPRODUCTION_TARGETS.csv')},indent=2))
            return int(args.require_models and (len(primary)!=3 or not all(x['available'] and assets.sha(root/x['path'])==x['sha256'] for x in primary)))
        if args.command=='smoke':
            if args.dry_run: print('Validated: synthetic CPU fixture; no research inputs or original states required');return 0
            from .workflow import smoke
            print(json.dumps(smoke(root,args.output),indent=2));return 0
        if args.command=='infer':
            from .workflow import infer
            print(json.dumps(infer(root,args.config,args.dry_run),indent=2));return 0
        if args.command in ['publication','build-documents']:
            from .publication import run,build
            return run(root,args) if args.command=='publication' else build(root,args)
    except (ValueError,FileNotFoundError,FileExistsError,KeyError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 2

if __name__=='__main__':sys.exit(main())
