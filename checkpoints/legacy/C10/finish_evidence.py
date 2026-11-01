"""Descriptive gene effects and package provenance; no model fitting/selection."""
from synthesis import *
import importlib.metadata

def main():
 g=pd.read_csv(B/'C10B_GENE_PATIENT_METRICS.csv');keys=['patient_id','regime','strategy','budget','gene'];wide=g.pivot(index=keys,columns='method',values=['pearson','mse','mae']);out=[]
 for a,b,comp in [('M1','M0','sparse_ST'),('M3','M1','image')]:
  d=wide.index.to_frame(index=False);d['comparison']=comp
  for m,sign in [('pearson',1),('mse',-1),('mae',-1)]:d[{'pearson':'delta_pearson','mse':'mse_reduction','mae':'mae_reduction'}[m]]=sign*(wide[m][a]-wide[m][b]).values
  out.append(d)
 eff=pd.concat(out,ignore_index=True);eff.to_csv(B/'C10B_GENE_PATIENT_GAINS.csv',index=False)
 q=eff[(eff.regime=='V1')&(eff.strategy=='random')&(eff.budget==.2)];rows=[]
 for key,x in q.groupby(['gene','comparison']):rows.append({'gene':key[0],'comparison':key[1],'median_delta_pearson':x.delta_pearson.median(),'median_mse_reduction':x.mse_reduction.median(),'median_mae_reduction':x.mae_reduction.median(),'patients_Pearson_improved':int((x.delta_pearson>0).sum()),'patients_all_metrics_improved':int(((x.delta_pearson>0)&(x.mse_reduction>0)&(x.mae_reduction>0)).sum()),'patients_defined':int(x.delta_pearson.notna().sum()),'source_file':'C10/C10B/C10B_GENE_PATIENT_GAINS.csv','source_rows':','.join(map(str,x.index+2))})
 pd.DataFrame(rows).to_csv(B/'C10B_GENE_GAIN_SUMMARY.csv',index=False)
 env=json.loads((B/'logs/ENVIRONMENT.json').read_text());env['versions']={name:importlib.metadata.version(name) for name in ['torch','torchvision','numpy','scipy','pandas','matplotlib','h5py','pillow','threadpoolctl']};dump(B/'logs/ENVIRONMENT.json',env)
 # Add clean topology controls to diagnostic source; no old experiment rerun.
 sources=json.loads((M/'SOURCE_INDEX.json').read_text());d=read('C3/C3A/C3A_PANEL_METRICS.csv');save(d,'SFig3_panel_source');d=read('C10/C10A_DATASET_MANIFEST.csv');save(d,'SFig10_QC_source');sources.update(SOURCES);dump(M/'SOURCE_INDEX.json',sources)
 print('GENE EFFECTS',len(eff),len(rows))
if __name__=='__main__':main()
