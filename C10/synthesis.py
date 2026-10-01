"""Read-only synthesis of frozen checkpoints. No prediction/training changes."""
from contract10 import *
M=OUT/'MANUSCRIPT'
for name in ['figures','supp_figures','tables','supp_tables','data']:(M/name).mkdir(parents=True,exist_ok=True)
SOURCES={}
def read(name):
 p=ROOT/name;d=pd.read_csv(p);d['source_file']=name;d['source_row']=np.arange(len(d))+2;d['source_sha256']=sha(p);SOURCES[name]={'sha256':sha(p),'rows':len(d)};return d

def agg(d,keys,metrics):
 rows=[]
 for key,g in d.groupby(keys,dropna=False,sort=True):
  if len(keys)==1:key=(key,)
  rows.append({**dict(zip(keys,key)),**{m:g[m].mean() for m in metrics},'source_file':'|'.join(g.source_file.unique()),'source_rows':','.join(map(str,g.source_row))})
 return pd.DataFrame(rows)

def save(d,name,folder='data',markdown=False):
 p=M/folder/(name+'.csv');d.to_csv(p,index=False)
 if markdown:(p.with_suffix('.md')).write_text('# '+name.replace('_',' ')+'\n\n'+('Full machine-readable table: '+p.name+'. Markdown shows the first 100 rows.\n\n' if len(d)>100 else '')+d.head(100).to_markdown(index=False)+'\n')
 return d

def main():
 rows=[];metrics=['median_pearson','mean_mse','mean_mae']
 c5=read('C5/C5_PANEL_METRICS.csv');sel=json.loads((ROOT/'C5/SELECTION_LOCK.json').read_text())['simple']['S2-50'];c5=c5[(c5.phase=='primary')&(c5.panel=='S2-50')&c5.control.isna()]
 c5=c5[((c5.model=='B1')|(c5.model=='B3'))&(c5.encoder==ENC)|((c5.encoder=='none')&c5.apply(lambda r:r.model==sel[r.regime]['best_family'],axis=1))].copy()
 c5['method']=c5.model.map({'B1':'M0','B2a':'M1','B2b':'M1','B3':'M3'});c5['regime']=c5.regime.map({'random':'V1','block':'V2','capture_window':'capture_window'});c5=c5.rename(columns={'nominal_coverage':'coverage','sample_id':'patient'})
 dev=agg(c5,['patient','coverage','regime','method'],metrics);save(c5,'Fig2_source');save(dev,'development_patient')
 def add(d,cohort,study,role,platform,ng,notes):
  for _,r in d.iterrows():rows.append({**r.to_dict(),'cohort':cohort,'study':study,'role':role,'platform':platform,'primary_gene_count':ng,'distance_metric_summary':'See Supplementary Table 7; occupied-bin descriptive curves, no fixed horizon','notes':notes})
 add(dev,'Development','C0 selected breast Xenium studies','development_primary','Xenium',50,'One primary patient; mean of3masks. V2 C5block differs from external V2. C5 fixed fusion weights differ from C9/C10.')
 hist=read('C9/C9_TRANSPORTABILITY_SUMMARY.csv');her=hist[hist.cohort.str.startswith('HER2ST')].copy().rename(columns={'mask_regime':'regime','mse':'mean_mse','mae':'mean_mae'});save(her,'Fig4_source')
 for old,new,methods,role in [('HER2ST_C8A','HER2ST zero-shot',{'T0_E0':'M0','T0_E1':'M1','T0_E2':'M3'},'transportability_development'),('HER2ST_C8A','HER2ST calibrated',{'T1_E0':'M0','T1_E1':'M1','T1_E2':'M3'},'transportability_development'),('HER2ST_C8B','HER2ST selected fusion',{'T1_F0':'M1','T1_SELECTED':'M3'},'transportability_development')]:
  d=her[(her.cohort==old)&her.method.isin(methods)].copy();d['method']=d.method.map(methods);d['source_rows']=d.source_row.astype(str);add(d[['patient','coverage','regime','method']+metrics+['source_file','source_rows']],new,'Andersson et al. HER2ST','transportability_development','Spatial Transcriptomics array',40,'Exact common40genes; not fresh validation. T0 replays C7. T1 errors in common C0 units. Affine interpolation equivalent.')
 for checkpoint,cohort,study in [('C9','Wu fresh validation','Wu et al. 2021'),('C10','Bassiouni second replication','Bassiouni et al. 2023; GSE210616')]:
  d=read(f'{checkpoint}/{checkpoint}B/{checkpoint}B_PATIENT_METRICS.csv');d=d[(d.strategy=='random')&d.method.isin(['M0','M1','M3'])].copy().rename(columns={'patient_id':'patient','budget':'coverage'});d['source_rows']=d.source_row.astype(str);save(d,'Fig5_source' if checkpoint=='C9' else 'Fig6_source')
  add(d[['patient','coverage','regime','method']+metrics+['source_file','source_rows']],cohort,study,'fresh_external_validation','Visium',50,'Patient mean over10seeds and equal V2geometries; median acrosspatients. Public H&E downsampled. '+('3heldout of6;3calibration' if checkpoint=='C9' else '22heldout, no C10calibration; M3 transports frozen Wu scaler.'))
 cross=pd.DataFrame(rows);cols=['cohort','study','patient','role','platform','primary_gene_count','coverage','regime','method']+metrics+['distance_metric_summary','notes','source_file','source_rows'];cross=cross[cols];cross.to_csv(OUT/'CROSS_COHORT_SUMMARY.csv',index=False);save(cross,'SuppTable2_per_patient','supp_tables',True)
 # Every descriptive aggregate retains explicit upstream CSV row numbers.
 tables=[]
 for key,g in cross.groupby(['cohort','primary_gene_count','coverage','regime','method']):
  row=dict(zip(['cohort','primary_gene_count','coverage','regime','method'],key));support=0
  for _,r in g.iterrows():
   base=cross[(cross.cohort==r.cohort)&(cross.patient==r.patient)&(cross.coverage==r.coverage)&(cross.regime==r.regime)&(cross.method=='M0')]
   if len(base)==1:support+=int(r.median_pearson>base.iloc[0].median_pearson and r.mean_mse<base.iloc[0].mean_mse and r.mean_mae<base.iloc[0].mean_mae)
  tables.append({**row,**{m:g[m].median() for m in metrics},'patient_count':g.patient.nunique(),'patient_support_vs_M0_all_metrics':support,'aggregation':'Median of patient summaries; masks/seeds averaged within patient','source_file':'C10/CROSS_COHORT_SUMMARY.csv','source_rows':','.join(map(str,g.index+2))})
 table2=save(pd.DataFrame(tables),'Table2_reconstruction','tables',True)
 t3=[]
 for key,g in cross[(cross.method=='M1')&cross.coverage.isin([.1,.2,.4])].groupby(['cohort','primary_gene_count']):
  r={'cohort':key[0],'primary_gene_count':key[1]}
  for b in [.1,.2,.4]:
   x=g[(g.regime=='V1')&(g.coverage==b)]
   for met in metrics[:2]:r[f'V1_{int(b*100)}pct_{met}']=x[met].median()
  x=g[g.coverage==.2].pivot(index='patient',columns='regime',values=metrics)
  if ('median_pearson','V2') in x:
   for met in metrics[:2]:r['V2_20pct_'+met]=x[met]['V2'].median();r['paired_V2_minus_V1_'+met]=(x[met]['V2']-x[met]['V1']).median()
  r['source_file']='C10/CROSS_COHORT_SUMMARY.csv';r['source_rows']=','.join(map(str,g.index+2));t3.append(r)
 save(pd.DataFrame(t3),'Table3_coverage_geometry','tables',True)
 # Distance data exported exactly, including location counts and sparse bins.
 ds=[]
 d=read('C5/C5_DISTANCE_TO_ANCHOR.csv');d=d[d.run_id.isin(c5.run_id)];d['cohort']='Development';d['method']=d.model.map({'B1':'M0','B2a':'M1','B2b':'M1','B3':'M3'});d['regime']=d.regime.map({'random':'V1','block':'V2','capture_window':'capture_window'});d=d.rename(columns={'nominal_coverage':'coverage','sample_id':'patient','replicate':'geometry'});ds.append(d)
 for cp,label in [('C9','Wu fresh validation'),('C10','Bassiouni second replication')]:
  d=read(f'{cp}/{cp}B/{cp}B_DISTANCE_TO_ANCHOR.csv');d=d[d.control.isna()&(d.strategy=='random')&d.method.isin(['M0','M1','M3'])].copy();d['cohort']=label;d=d.rename(columns={'budget':'coverage','patient_id':'patient'});ds.append(d)
 distances=pd.concat(ds,ignore_index=True);save(distances,'SuppTable7_distance_bins','supp_tables',True)
 # per-gene fresh replication metrics plus original historical tables, unmodified values.
 parts=[]
 for cp in ['C9','C10']:
  d=read(f'{cp}/{cp}B/{cp}B_GENE_METRICS.csv.gz' if cp=='C10' else f'{cp}/{cp}B/{cp}B_GENE_METRICS.csv');d['checkpoint']=cp;parts.append(d)
 gene=pd.concat(parts,ignore_index=True);gene.to_csv(M/'supp_tables/SuppTable3_per_gene.csv.gz',index=False,compression='gzip');(M/'supp_tables/SuppTable3_per_gene.md').write_text('# Per-gene metrics\n\nComplete fresh C9/C10 metric rows are in `SuppTable3_per_gene.csv.gz`; no rounding. Historical C1–C8 gene metrics remain frozen at paths in SOURCE_INDEX.json. Patient/seed/geometry/gene are identified; these rows are not independent biological replicates.\n\n'+gene.head(25).to_markdown(index=False)+'\n')
 g=gene[(gene.checkpoint=='C10')&gene.control.isna()];g.groupby(['patient_id','regime','strategy','budget','method','gene'])[['pearson','spearman','mse','mae']].mean().reset_index().to_csv(B/'C10B_GENE_PATIENT_METRICS.csv',index=False)
 # Exact panels including missing symbols.
 panels=[]
 for panel,gg in parent.ORIGINAL.items():
  for symbol in gg:panels.append({'checkpoint':'C0','panel':panel,'gene':symbol,'included':True})
 for cp,path in [('C10','C10/C10A/GENE_INTERSECTION.csv')]:
  d=read(path);d['checkpoint']=cp;panels.extend(d.to_dict('records'))
 d=read('C9/C9A/GENE_INTERSECTION.csv') if (ROOT/'C9/C9A/GENE_INTERSECTION.csv').exists() else None
 if d is not None:d['checkpoint']='C9';panels.extend(d.to_dict('records'))
 save(pd.DataFrame(panels),'SuppTable1_gene_panels','supp_tables',True)
 # Supplementary plotting sources retain exact frozen source rows.
 names=[('SFig1','C1/C1_METRICS_SUMMARY.csv'),('SFig1_controls','C1/C1_NEGATIVE_CONTROLS.csv'),('SFig2','C2/C2_METRICS_SUMMARY.csv'),('SFig2_controls','C2/C2_GRAPH_CONTROLS.csv'),('SFig3','C3/C3A/C3A_RANDOMIZATION_METRICS.csv'),('SFig4','C4/C4A/C4A_PANEL_METRICS.csv'),('SFig5','C5/C5_NEGATIVE_CONTROLS.csv'),('SFig6','C6/C6A/C6A_PANEL_METRICS.csv'),('SFig7','C7/C7B/C7B_COVERAGE_CURVES.csv'),('SFig8','C8/C8B/C8B_PATIENT_EFFECTS.csv'),('SFig9','C9/C9B/C9B_NEGATIVE_CONTROLS.csv'),('SFig10','C10/C10B/C10B_NEGATIVE_CONTROLS.csv')]
 nc=[]
 for name,path in names:
  d=read(path);save(d,name+'_source');
  if 'control' in name or name in ['SFig5','SFig9','SFig10']:d['source_checkpoint']=path.split('/')[0];nc.append(d)
 save(pd.concat(nc,ignore_index=True),'SuppTable4_negative_controls','supp_tables',True)
 # Source index also points to every frozen historical per-gene table.
 for cp in range(1,9):
  for p in (ROOT/f'C{cp}').rglob('*GENE*csv*'):
   if p.name.endswith('.csv') and p.with_suffix('.csv.gz').exists():continue
   SOURCES[str(p.relative_to(ROOT))]={'sha256':sha(p),'use':'Historical per-gene evidence, unchanged; see checkpoint reports for estimand'}
 dump(M/'SOURCE_INDEX.json',SOURCES)
 print('SYNTHESIS',len(cross),'patient rows',len(table2),'table rows',flush=True)
if __name__=='__main__':main()
