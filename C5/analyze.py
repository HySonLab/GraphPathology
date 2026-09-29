"""All hidden-location metrics and predeclared sparse-reconstruction decision rules."""
from common import *
import math,gzip
from scipy.stats import spearmanr

def main():
 records=[json.loads(p.read_text()) for p in sorted((OUT/'logs/runs').glob('*.json'))];assert len(records)==3456,len(records)
 table=pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in records]);table.to_csv(OUT/'C5_PANEL_METRICS.csv',index=False);table[table.control.notna()].to_csv(OUT/'C5_NEGATIVE_CONTROLS.csv',index=False);table[table.phase=='sensitivity'].to_csv(OUT/'C5_PATIENT_SENSITIVITY.csv',index=False)
 arrays={};gene_frames=[];dist=[]
 for r in records:
  with np.load(ROOT/r['prediction_path']) as f:m={k:f[k].copy() for k in ['pearson','spearman','mse','mae']};genes=f['genes'].copy()
  arrays[r['run_id']]=m;cols={k:r[k] for k in ['run_id','cohort','sample_id','panel','mask_id','regime','nominal_coverage','realized_coverage','replicate','model','encoder','control','seed']};gene_frames.append(pd.DataFrame({'gene':genes,**m}).assign(**cols))
  dist.extend([{**cols,**d} for d in r['distance_bins']])
 genes=pd.concat(gene_frames,ignore_index=True);genes.to_csv(OUT/'C5_GENE_METRICS.csv',index=False)
 with (OUT/'C5_GENE_METRICS.csv').open('rb') as src,gzip.open(OUT/'C5_GENE_METRICS.csv.gz','wb') as dst:
  import shutil;shutil.copyfileobj(src,dst)
 pd.DataFrame(dist).to_csv(OUT/'C5_DISTANCE_TO_ANCHOR.csv',index=False)
 ordinary=table[table.control.isna()];group=['cohort','sample_id','panel','regime','nominal_coverage','model','encoder'];curves=ordinary.groupby(group).agg(realized_coverage=('realized_coverage','mean'),median_pearson=('median_pearson','mean'),pearson_mask_SD=('median_pearson','std'),mean_mse=('mean_mse','mean'),mse_mask_SD=('mean_mse','std'),mean_mae=('mean_mae','mean'),n_masks=('run_id','size')).reset_index();curves.to_csv(OUT/'C5_COVERAGE_CURVES.csv',index=False)
 gc=genes[genes.control.isna()].groupby(group+['gene'])[['realized_coverage','pearson','spearman','mse','mae']].mean().reset_index();gc.to_csv(OUT/'C5_GENE_COVERAGE.csv',index=False)
 with (OUT/'C5_GENE_COVERAGE.csv').open('rb') as src,gzip.open(OUT/'C5_GENE_COVERAGE.csv.gz','wb') as dst:shutil.copyfileobj(src,dst)
 selection=json.loads((OUT/'SELECTION_LOCK.json').read_text());lookup={(r['cohort'],r['panel'],r['mask_id'],r['model'],r['encoder']):r for r in records if r['control'] is None};settings=[];gains=[]
 for key,d in ordinary[ordinary.model=='B3'].groupby(['cohort','sample_id','panel','regime','nominal_coverage','encoder']):
  cohort,sample,panel,regime,coverage,enc=key;sp=[];fusion=[];mse=[];b4gain=[];b4mse=[];ability=[]
  for r in d.sort_values('replicate').to_dict('records'):
   get=lambda model,encoder:lookup[(cohort,panel,r['mask_id'],model,encoder)];b2=get(selection['simple'][panel][regime]['best_family'],'none');b1=get('B1',enc);b0=get('B0','none');b4=get('B4',enc);a2=arrays[b2['run_id']];a3=arrays[r['run_id']];a4=arrays[b4['run_id']]
   sp.append(b2['median_pearson']-b1['median_pearson']>=.02 and b2['mean_mse']<=.98*b1['mean_mse']);fusion.append(a3['pearson']-a2['pearson']);mse.append(a2['mse']-a3['mse']);b4gain.append(a4['pearson']-a3['pearson']);b4mse.append(a3['mse']-a4['mse']);ability.append(r['median_pearson']>=.2 and r['mean_mse']<=.9*b0['mean_mse'])
   gains.append(pd.DataFrame({'gene':PANELS[panel],'H_E_Pearson_gain':fusion[-1],'H_E_MSE_reduction':mse[-1],'propagation_Pearson_gain':b4gain[-1],'propagation_MSE_reduction':b4mse[-1]}).assign(**{k:r[k] for k in ['cohort','sample_id','panel','mask_id','regime','nominal_coverage','encoder']}))
  f=np.array(fusion);e=np.array(mse);b=np.array(b4gain);be=np.array(b4mse);strong=int(((np.nanmean(f,0)>=.05)&((f>0).sum(0)>=2)&(e.mean(0)>=0)).sum());qualified=bool(np.nanmedian(np.nanmean(f,0))>=.01 and e.mean()>=0 and strong>=max(3,math.ceil(.05*f.shape[1])))
  settings.append({'cohort':cohort,'sample_id':sample,'panel':panel,'regime':regime,'nominal_coverage':coverage,'encoder':enc,'sparse_setting':sum(sp)>=2,'fusion_setting':qualified,'strong_fusion_genes':strong,'median_fusion_gene_gain':float(np.nanmedian(np.nanmean(f,0))),'mean_fusion_MSE_reduction':float(e.mean()),'contiguous_ability':regime!='random' and sum(ability)>=2,'propagation_setting':bool(np.nanmedian(np.nanmean(b,0))>=.01 and be.mean()>=0 and (np.nanmedian(b,1)>0).sum()>=2),'median_propagation_gene_gain':float(np.nanmedian(np.nanmean(b,0))),'mean_propagation_MSE_reduction':float(be.mean())})
 settings=pd.DataFrame(settings);settings.to_csv(OUT/'C5_SETTING_GAINS.csv',index=False);pd.concat(gains,ignore_index=True).to_csv(OUT/'C5_GENE_GAINS.csv.gz',index=False,compression='gzip')
 coherence=[]
 for panel in PANELS:
  q=curves[(curves.cohort=='primary')&(curves.panel==panel)&(curves.regime=='random')];model=selection['simple'][panel]['random']['best_family'];q=q[q.model==model].sort_values('realized_coverage');rp=float(spearmanr(q.realized_coverage,q.median_pearson).statistic);re=float(spearmanr(q.realized_coverage,q.mean_mse).statistic);coherence.append({'panel':panel,'Pearson_coverage_rho':rp,'MSE_coverage_rho':re,'coherent':rp>0 and re<0,'MSE_low_to_high_relative_reduction':float(1-q.mean_mse.iloc[-1]/q.mean_mse.iloc[0])})
 rows=[]
 for panel in PANELS:
  for enc in ENCODERS:
   p=settings[(settings.cohort=='primary')&(settings.panel==panel)&(settings.encoder==enc)];f=settings[(settings.cohort!='primary')&(settings.panel==panel)&(settings.encoder==enc)];patients=0;contiguous=False
   for regime in REGIMES:
    primary=p[p.regime==regime];fold=f[f.regime==regime];count=sum(bool(d.fusion_setting.any() and d.sparse_setting.any()) for _,d in fold.groupby('sample_id'))
    if primary.fusion_setting.any():patients=max(patients,count)
    if regime!='random' and primary.contiguous_ability.any() and sum(d.contiguous_ability.any() for _,d in fold.groupby('sample_id'))>=2:contiguous=True
   passed=bool(p.sparse_setting.sum()>=2 and p.fusion_setting.sum()>=2 and patients>=2 and contiguous and sum(r['coherent'] for r in coherence)>=2)
   rows.append({'panel':panel,'encoder':enc,'primary_sparse_settings':int(p.sparse_setting.sum()),'primary_fusion_settings':int(p.fusion_setting.sum()),'supporting_patients_same_geometry':int(patients),'contiguous_ability_supported':bool(contiguous),'pass':passed})
 status='PASS_SPARSE_ST_RECONSTRUCTION' if any(r['pass'] for r in rows) else ('MIXED_SPARSE_ST_RECONSTRUCTION' if settings[settings.cohort=='primary'][['sparse_setting','fusion_setting']].any().any() else 'FAIL_SPARSE_ST_RECONSTRUCTION')
 recommendation=('PROCEED_TO_SAMPLING_OPTIMIZATION' if any(r['MSE_low_to_high_relative_reduction']>=.1 for r in coherence) else 'PROCEED_TO_EXTERNAL_VALIDATION') if status=='PASS_SPARSE_ST_RECONSTRUCTION' else 'NEED_SPARSE_ST_REVIEW'
 dump(OUT/'logs/GATE_RESULT.json',{'status':status,'recommendation':recommendation,'panel_encoder_rules':rows,'coherence':coherence,'technical_audit_required':True,'prediction_runs':len(records),'settings':len(settings)})
 events=[json.loads(f.read_text()) for f in (OUT/'logs/registry').glob('*.json')];pd.DataFrame([{k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in events]).to_csv(OUT/'C5_REGISTRY_EXTENSIONS.csv',index=False);table[['run_id','cohort','panel','model','encoder','n_hidden','inference_audit_metric_seconds','peak_rss_mb']].to_csv(OUT/'C5_RESOURCE_PROFILE.csv',index=False);print(status,recommendation,flush=True)
if __name__=='__main__':main()
