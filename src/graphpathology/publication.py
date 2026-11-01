"""Regenerate available export views, with explicit revised-source limitations."""
import csv,json,shutil
from pathlib import Path
from .assets import new_output,sha

def targets(root):return list(csv.DictReader((root/'docs/REPRODUCTION_TARGETS.csv').open()))
def source(root,name):return root/'checkpoints/legacy'/name

def run(root,args):
    rows=targets(root)
    if args.target!='all':rows=[x for x in rows if x['target']==args.target]
    if not rows:raise ValueError('Unknown target; consult docs/REPRODUCTION_TARGETS.csv')
    blocked=[x for x in rows if x['status'].startswith('blocked')]
    if args.dry_run:
        print(json.dumps({'targets':rows,'revised_baseline_verified':False},indent=2));return 2 if blocked else 0
    if blocked and not args.available_only:raise FileNotFoundError('Revised inputs missing for '+', '.join(x['target'] for x in blocked)+'; obtain revised LaTeX/scripts/source_data/REPRODUCTION_TARGETS.csv package from author. Use --available-only for labeled older-export views.')
    if not args.output:raise ValueError('--output is required')
    out=new_output(root,args.output)
    import numpy as np,pandas as pd
    from . import export_plots as plot
    plot.ROOT=root/'checkpoints/legacy';plot.M=out;plot.D=out/'data'
    for name in ['data','figures','supp_figures','tables','supp_tables']:(out/name).mkdir()
    baseline=source(root,'C10/MANUSCRIPT')
    # Input bytes are preserved in the new work tree; never mutate the archived package.
    for name in ['tables','supp_tables']:
        for p in (baseline/name).glob('*.csv*'):shutil.copy2(p,out/name/p.name)
    def save(frame,name):frame.to_csv(out/'data'/(name+'.csv'),index=False)
    cross=pd.read_csv(source(root,'C10/CROSS_COHORT_SUMMARY.csv'))
    save(cross[cross.cohort=='Development'],'development_patient')
    save(cross[cross.cohort=='Wu fresh validation'],'Fig5_source')
    save(cross[cross.cohort=='Bassiouni second replication'],'Fig6_source')
    her=pd.read_csv(source(root,'C9/C9_TRANSPORTABILITY_SUMMARY.csv')).rename(columns={'mask_regime':'regime','mse':'mean_mse','mae':'mean_mae'});save(her,'Fig4_source')
    names={'SFig1':'C1/C1_METRICS_SUMMARY.csv','SFig1_controls':'C1/C1_NEGATIVE_CONTROLS.csv','SFig2':'C2/C2_METRICS_SUMMARY.csv','SFig2_controls':'C2/C2_GRAPH_CONTROLS.csv','SFig3':'C3/C3A/C3A_RANDOMIZATION_METRICS.csv','SFig3_panel':'C3/C3A/C3A_PANEL_METRICS.csv','SFig4':'C4/C4A/C4A_PANEL_METRICS.csv','SFig5':'C5/C5_NEGATIVE_CONTROLS.csv','SFig6':'C6/C6A/C6A_PANEL_METRICS.csv','SFig7':'C7/C7B/C7B_COVERAGE_CURVES.csv','SFig8':'C8/C8B/C8B_PATIENT_EFFECTS.csv','SFig9':'C9/C9B/C9B_NEGATIVE_CONTROLS.csv','SFig10':'C10/C10B/C10B_NEGATIVE_CONTROLS.csv','SFig10_QC':'C10/C10A_DATASET_MANIFEST.csv'}
    for name,p in names.items():save(pd.read_csv(source(root,p)),name+'_source')
    # S7 was revised to a single random mask, rather than the old all-geometry plot.
    p=out/'data/SFig5_source.csv';d=pd.read_csv(p);d=d[(d.regime=='random')&(d.replicate==0)];d.to_csv(p,index=False)
    results=[]
    for row in rows:
        target=row['target']
        if row in blocked:results.append({'target':target,'status':row['status']});continue
        before=set(out.rglob('*'))
        if target=='Figure3':plot.fig2()
        elif target=='Figure4':plot.external(5,'Wu fresh validation')
        elif target=='Figure5':plot.external(6,'Bassiouni second replication')
        elif target=='FigureS1':plot.fig4()
        elif target=='FigureS2':plot.fig7()
        elif target.startswith('FigureS'):plot.supplement(int(target[7:])-2)
        else:
            if target=='TableS8':
                from .graph import parameter_count
                records=json.loads((root/'manuscript/tables/TableS8.json').read_text())['rows']
                for rec in records:
                    expected=parameter_count(rec['input'],rec['hidden'],rec['output']) if rec['hidden'] is not None else (rec['input']+1)*rec['output']
                    if expected!=rec['parameters']:raise ValueError('TableS8 parameter count mismatch: '+rec['predictor'])
                table=pd.DataFrame(records)
            elif target=='TableS9':
                d=cross[(cross.method=='M1')&cross.cohort.isin(['Wu fresh validation','Bassiouni second replication'])&cross.coverage.isin([.1,.2,.4])]
                table=d.groupby(['cohort','coverage','regime'])[['median_pearson','mean_mse']].median().reset_index()
            else:
                p=source(root,row['source']);table=pd.read_csv(p,low_memory=False)
                if target=='Table2':table=table[(table.coverage==.2)&(table.regime=='V1')&table.cohort.isin(['Development','HER2ST zero-shot','Wu fresh validation','Bassiouni second replication'])&table.method.isin(['M0','M1','M3'])]
                if target=='TableS1':
                    genes=json.loads(source(root,'reports/C0/C0_GENE_PANELS.json').read_text())
                    # Ordered primary gene list is canonical in the external protocol.
                    genes=json.loads(source(root,'C9/configs/protocol.json').read_text())['genes']
                    table=pd.DataFrame({'rank':np.arange(1,51),'gene':genes})
                if target=='TableS2':table=table[(table.coverage==.2)&(table.regime=='V1')&table.cohort.isin(['Wu fresh validation','Bassiouni second replication'])]
            table.to_csv(out/'tables'/(target+'.csv'),index=False)
            # Large gene/mask exports remain complete CSVs; previews are explicitly limited.
            (out/'tables'/(target+'.html')).write_text('<meta charset="utf-8"><h1>'+target+'</h1><p>Export view; revised TeX template unavailable. First 100 rows shown; complete CSV retained.</p>'+table.head(100).to_html(index=False))
        created=[p for p in set(out.rglob('*'))-before if p.is_file()]
        # Rename figure files to revised identifiers; these are export views, not a revised document build.
        for p in created:
            if p.parent.name in ['figures','supp_figures']:
                dest=p.parent/(target+p.suffix);p.rename(dest)
            elif p.name.endswith('_plotted_values.json'):p.rename(out/(target+'_plotted_values.json'))
        results.append({'target':target,'status':'export_view_generated','revised_baseline_verified':False})
    traces=[]
    for p in sorted((out/'data').glob('*.csv')):
        d=pd.read_csv(p,low_memory=False);traces.append({'name':p.name,'sha256':sha(p),'rows':len(d),'columns':list(d.columns),'missing':{k:int(v) for k,v in d.isna().sum().items()}})
    (out/'VALUE_TRACES.json').write_text(json.dumps(traces,indent=2)+'\n')
    report={'targets':results,'revised_baseline_verified':False,'limitations':'Revised source package unavailable. These use older numerical exports and maintained renderers; no claim of full revised display equivalence or document compilation.'}
    (out/'TARGET_STATUS.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    return 0

def build(root,args):
    expected=root/'manuscript/revised/build.sh'
    if not expected.is_file():raise FileNotFoundError('Revised LaTeX source package missing: install under manuscript/revised/ (build.sh, TeX, bibliography, figures, source_data). Supplied PDFs cannot be compiled back into their sources.')
    raise ValueError('Revised source package newly installed: inspect its build script and update target manifest before executing; this checkout has not validated that build.')
