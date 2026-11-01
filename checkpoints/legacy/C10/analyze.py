"""Predeclared patient-level summaries and second-replication gates."""
from contract10 import *
import gzip

def main():
 check_lock();records=[json.loads(p.read_text()) for p in sorted((B/'logs/runs').glob('*.json'))];assert len(records)==len(SAMPLES)*132+18
 panels=[];dd=[];field=['run_id','patient_id','regime','geometry','strategy','seed','budget','control','control_seed'];gp=B/'C10B_GENE_METRICS.csv'
 with gp.open('w',newline='') as handle:
  writer=csv.DictWriter(handle,fieldnames=field+['method','gene','pearson','spearman','mse','mae']);writer.writeheader()
  for r in records:
   base={k:r[k] for k in field}
   with np.load(ROOT/r['prediction_path']) as f:
    for m in ['M0','M1','M3']:
     panels.append({**base,'method':m,**r['summary'][m]})
     for j,g in enumerate(genes()):writer.writerow({**base,'method':m,'gene':g,**{k:float(f[m+'_'+k][j]) for k in ['pearson','spearman','mse','mae']}})
   dd.extend({**base,**d} for d in r['distance_bins'])
 with gp.open('rb') as src,gzip.open(gp.with_suffix('.csv.gz'),'wb') as dst:
  import shutil;shutil.copyfileobj(src,dst)
 df=pd.DataFrame(panels);df.to_csv(B/'C10B_PANEL_METRICS.csv',index=False);pd.DataFrame(dd).to_csv(B/'C10B_DISTANCE_TO_ANCHOR.csv',index=False);clean=df[df.control==''];metrics=list(records[0]['summary']['M0']);keys=['patient_id','regime','strategy','budget','method'];geo=clean.groupby(keys+['geometry'])[metrics].mean().reset_index();pt=geo.groupby(keys)[metrics].mean().reset_index();pt.to_csv(B/'C10B_PATIENT_METRICS.csv',index=False);geo.to_csv(B/'C10B_GEOMETRY_METRICS.csv',index=False);pt.to_csv(B/'C10B_COVERAGE_CURVES.csv',index=False)
 effects=[]
 for key,g in pt.groupby(['patient_id','regime','strategy','budget']):
  g=g.set_index('method');base=dict(zip(['patient_id','regime','strategy','budget'],key))
  for a,b,label in [('M1','M0','sparse_ST'),('M3','M1','image')]:
   dr=float(g.loc[a,'median_pearson']-g.loc[b,'median_pearson']);mse=float(g.loc[b,'mean_mse']-g.loc[a,'mean_mse']);mae=float(g.loc[b,'mean_mae']-g.loc[a,'mean_mae']);rel=mse/g.loc[b,'mean_mse']
   effects.append({**base,'comparison':label,'delta_pearson':dr,'mse_reduction':mse,'relative_mse_reduction':rel,'mae_reduction':mae,'meaningful_sparse':dr>=.02 and rel>=.02 and mae>=0})
 ef=pd.DataFrame(effects);ef.to_csv(B/'C10B_PATIENT_EFFECTS.csv',index=False)
 nc=df[df.control!=''].copy();join=['patient_id','regime','geometry','strategy','seed','budget','method'];matched=clean[join+metrics].rename(columns={k:'clean_'+k for k in metrics});nc=nc.merge(matched,on=join,validate='many_to_one');nc['delta_pearson']=nc.median_pearson-nc.clean_median_pearson;nc['mse_reduction']=nc.clean_mean_mse-nc.mean_mse;nc['mae_reduction']=nc.clean_mean_mae-nc.mean_mae;nc.to_csv(B/'C10B_NEGATIVE_CONTROLS.csv',index=False)
 primary=ef[(ef.regime=='V1')&(ef.strategy=='random')];st=primary[primary.comparison=='sparse_ST'];support=st.groupby('patient_id').meaningful_sparse.sum();coverage=[];degradation=[]
 for s in SAMPLES:
  p=pt[(pt.patient_id==s)&(pt.strategy=='random')&(pt.method=='M1')];v1=p[p.regime=='V1'].set_index('budget');v2=p[p.regime=='V2'].set_index('budget')
  coverage.append({'patient_id':s,'delta_pearson_40_minus_10':float(v1.loc[.4,'median_pearson']-v1.loc[.1,'median_pearson']),'mse_reduction_10_to_40':float(v1.loc[.1,'mean_mse']-v1.loc[.4,'mean_mse'])})
  for b in [.1,.2,.4]:degradation.append({'patient_id':s,'budget':b,'delta_pearson_V2_minus_V1':float(v2.loc[b,'median_pearson']-v1.loc[b,'median_pearson']),'mse_increase_V2_minus_V1':float(v2.loc[b,'mean_mse']-v1.loc[b,'mean_mse']),'relative_mse_increase':float(v2.loc[b,'mean_mse']/v1.loc[b,'mean_mse']-1),'mae_increase':float(v2.loc[b,'mean_mae']-v1.loc[b,'mean_mae'])})
 cov=pd.DataFrame(coverage);cov.to_csv(B/'C10B_COVERAGE_EFFECTS.csv',index=False);pd.DataFrame(degradation).to_csv(B/'C10B_V2_DEGRADATION.csv',index=False);budget_gains=st.groupby('budget')[['delta_pearson','mse_reduction','mae_reduction']].median();required=int(np.ceil(2*len(SAMPLES)/3));coverage_patients=int(((cov.delta_pearson_40_minus_10>0)&(cov.mse_reduction_10_to_40>0)).sum());core_ok=bool((support>=2).sum()>=required and (budget_gains>0).all().all() and coverage_patients>=required)
 image=primary[primary.comparison=='image'];imok=(image.delta_pearson>0)&(image.mse_reduction>0)&(image.mae_reduction>0);image_ok=bool(imok.groupby(image.patient_id).sum().ge(2).sum()>=required and (image.relative_mse_reduction<-.02).mean()<=.25)
 status='PASS_LIMITED_SECOND_REPLICATION' if core_ok else ('MIXED_SECOND_REPLICATION' if int(st[st.meaningful_sparse].patient_id.nunique())>=2 else 'FAIL_SECOND_REPLICATION')
 dump(B/'logs/SCIENTIFIC_GATE.json',{'status_pending_technical_audit':status,'core_sparse_support':core_ok,'required_patients':required,'patient_meaningful_budget_counts':support.to_dict(),'coverage_response_patients':coverage_patients,'fusion_secondary_pass':image_ok,'runs':len(records),'method_predictions':len(records)*3,'corruption_invariance_checks':sum(r['hidden_corruption_checks'] for r in records),'recommendation':'READY_WITH_LIMITED_SECOND_REPLICATION' if core_ok else 'NEED_METHOD_REVIEW'})
 print(json.dumps(json.loads((B/'logs/SCIENTIFIC_GATE.json').read_text()),indent=2),flush=True)
if __name__=='__main__':main()
