"""Re-summarize saved C5 gene metrics on the identical C7 gene intersections."""
from support import *
def main():
 rows=[];old=pd.read_csv(ROOT/'C5/C5_PANEL_METRICS.csv');old=old[(old.phase=='primary')&old.panel.isin(panels())&old.regime.isin(['random','block'])&old.nominal_coverage.isin([.1,.2,.4])&old.control.isna()]
 for _,r in old.iterrows():
  if r.model in ['B1','B3','B4']:
   if r.encoder!=ENC:continue
   model={'B1':'E0','B3':'E2','B4':'E3'}[r.model]
  elif r.model==SEL['simple'][r.panel][r.regime]['best_family']:model='E1'
  else:continue
  with np.load(ROOT/r.prediction_path) as f:
   idx=[list(f['genes']).index(g) for g in panels()[r.panel]];metrics={k:f[k][idx] for k in ['pearson','spearman','mse','mae','prediction_variance','target_variance']};m=core.summarize(metrics)
  rows.append({'dataset':'development_C5_primary','patient':r.sample_id,'panel':r.panel,'budget':r.nominal_coverage,'mask_regime':'V1' if r.regime=='random' else 'V2','model':model,'strategy':'random' if r.regime=='random' else 'C5_block','replicate':r.replicate,'median_pearson':m['median_pearson'],'mse':m['mean_mse'],'mae':m['mean_mae'],'genes':len(idx),'source_prediction_path':r.prediction_path,'geometry_note':'C5 original coverage-dependent blocks; not the C7 fixed forbidden regions' if r.regime=='block' else 'C5 saved random masks'})
 dev=pd.DataFrame(rows);dev=dev.groupby(['dataset','patient','panel','budget','mask_regime','model','strategy','genes','geometry_note'],as_index=False)[['median_pearson','mse','mae']].mean();ext=pd.read_csv(B/'C7B_COVERAGE_CURVES.csv');ext=ext.rename(columns={'patient_id':'patient','nominal_budget':'budget','regime':'mask_regime','mean_mse':'mse','mean_mae':'mae'});ext['dataset']='external_HER2ST';ext['genes']=ext.panel.map({p:len(g) for p,g in panels().items()});ext['geometry_note']=ext.mask_regime.map({'V1':'C7 all unobserved tissue','V2':'C7 three fixed 40% unavailable regions'});cols=list(dev.columns);both=pd.concat([dev,ext[cols]],ignore_index=True)
 for _,inds in both.groupby(['dataset','patient','panel','budget','mask_regime','strategy']).groups.items():
  sub=both.loc[inds].set_index('model');where=[i for i in inds if both.loc[i,'model']=='E2'];assert len(where)==1;i=where[0]
  both.loc[i,'fusion_gain_vs_image']=sub.loc['E2','median_pearson']-sub.loc['E0','median_pearson'];both.loc[i,'fusion_gain_vs_interpolation']=sub.loc['E2','median_pearson']-sub.loc['E1','median_pearson'];both.loc[i,'fusion_MSE_reduction_vs_image']=sub.loc['E0','mse']-sub.loc['E2','mse'];both.loc[i,'fusion_MSE_reduction_vs_interpolation']=sub.loc['E1','mse']-sub.loc['E2','mse']
 both.to_csv(OUT/'C7_DEVELOPMENT_EXTERNAL_COMPARISON.csv',index=False);print('COMPARISON',len(both),'rows; same frozen gene subsets; no parent model rerun')
if __name__=='__main__':main()
