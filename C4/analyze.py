"""Apply the pre-test C4A representation gate without tuning thresholds."""
from common import *
import math
CFG=json.loads((A/'configs/protocol.json').read_text())

def analyze():
 evals=[json.loads(p.read_text()) for p in (A/'logs/evaluation').glob('*.json')];p=pd.DataFrame(evals);g=pd.concat([pd.read_csv(ROOT/r['gene_metrics_path']) for r in evals],ignore_index=True)
 p.to_csv(A/'C4A_PANEL_METRICS.csv',index=False);g.to_csv(A/'C4A_GENE_METRICS.csv',index=False);p[p.phase=='negative'].to_csv(A/'C4A_NEGATIVE_CONTROLS.csv',index=False)
 clean=g[g.phase.isin(['primary','lopo'])];key=['phase','sample_id','panel','model','seed','gene'];w=clean.pivot(index=key,columns='encoder',values=['pearson','mse','mae']);w.columns=[f'{m}_{e}' for m,e in w.columns];w=w.reset_index();w['delta_encoder']=w.pearson_pathology_encoder_1-w.pearson_resnet18_reference;w['delta_mse']=w.mse_resnet18_reference-w.mse_pathology_encoder_1;w['relative_mse_increase']=-w.delta_mse/w.mse_resnet18_reference
 w.to_csv(A/'C4A_ENCODER_GAIN.csv',index=False);w[w.phase=='lopo'].to_csv(A/'C4A_PATIENT_SENSITIVITY.csv',index=False)
 means=w.groupby(['phase','sample_id','panel','model','gene']).agg(delta_encoder=('delta_encoder','mean'),minimum_seed_gain=('delta_encoder','min'),delta_mse=('delta_mse','mean'),mse_resnet=('mse_resnet18_reference','mean'),mse_pathology=('mse_pathology_encoder_1','mean')).reset_index();means['strong_gene']=(means.delta_encoder>=CFG['candidate_rule']['primary_mean_Pearson_gain_min'])&(means.minimum_seed_gain>0)&(means.delta_mse>=0);means.to_csv(A/'C4A_GENE_SUPPORT.csv',index=False)
 rule=CFG['gate'];negative=[]
 for model in CFG['models']:
  true=p[(p.phase=='primary')&(p.encoder=='pathology_encoder_1')&(p.panel=='S2-50')&(p.model==model)].median_pearson.mean()
  for control in ['NC1','NC2']:
   d=p[(p.phase=='negative')&(p.model==model)&(p.control==control)];assert len(d)==3
   value=d.median_pearson.mean();passed=bool(np.isfinite(value) and value<=rule['negative_controls_max_median_Pearson'] and true-value>=rule['negative_controls_min_drop_vs_unpermuted']);negative.append({'model':model,'control':control,'mean_median_Pearson':value,'unpermuted':true,'drop':true-value,'pass':passed})
 panelrows=[]
 for panel in CFG['panels'][1:]:
  d=means[(means.phase=='primary')&(means.panel==panel)&(means.model=='D1')];ridge=means[(means.phase=='primary')&(means.panel==panel)&(means.model=='D0')];q=p[(p.phase=='primary')&(p.panel==panel)&(p.model=='D1')].groupby('encoder')[['median_pearson','mean_mse']].mean()
  good=[];foldstats=[]
  for sample in IDS:
   f=means[(means.phase=='lopo')&(means.panel==panel)&(means.model=='D1')&(means.sample_id==sample)];gain=f.delta_encoder.median();relative=f.mse_pathology.mean()/f.mse_resnet.mean()-1;good.append(gain>0 and relative<=.02);foldstats.append(relative)
  r={'panel':panel,'median_panel_correlation_increase':q.loc['pathology_encoder_1','median_pearson']-q.loc['resnet18_reference','median_pearson'],'median_gene_gain':d.delta_encoder.median(),'mean_gene_gain':d.delta_encoder.mean(),'fraction_Pearson_improved':float((d.delta_encoder>0).mean()),'fraction_both_improved':float(((d.delta_encoder>0)&(d.delta_mse>0)).mean()),'primary_relative_MSE_increase':q.loc['pathology_encoder_1','mean_mse']/q.loc['resnet18_reference','mean_mse']-1,'strong_genes':int(d.strong_gene.sum()),'required_strong_genes':max(3,math.ceil(.05*len(d))),'conditional_patients_improved_in_both':int(sum(good)),'conditional_mean_relative_MSE_increase':float(np.mean(foldstats)),'ridge_median_gene_gain':ridge.delta_encoder.median()}
  flags=[r['median_panel_correlation_increase']>=rule['minimum_panel_median_Pearson_increase'],r['median_gene_gain']>=rule['minimum_median_paired_gene_gain'],r['primary_relative_MSE_increase']<=rule['maximum_primary_mean_MSE_relative_increase'],r['fraction_both_improved']>=rule['minimum_fraction_genes_improved_in_both'],r['strong_genes']>=r['required_strong_genes'],r['conditional_patients_improved_in_both']>=rule['minimum_conditional_patients_with_positive_median_gain_and_MSE_increase_at_most_2percent'],r['conditional_mean_relative_MSE_increase']<=rule['maximum_conditional_average_MSE_relative_increase'],r['ridge_median_gene_gain']>=rule['same_panel_Ridge_median_gene_gain_at_least']]
  r['pass']=all(flags);panelrows.append(r)
 if not all(r['pass'] for r in negative):status='FAIL_TECHNICAL'
 elif any(r['pass'] for r in panelrows):status='PASS_PATHOLOGY_REPRESENTATION'
 elif any(r['median_gene_gain']>=.01 or r['strong_genes']>=.05*len(PANELS[r['panel']]) for r in panelrows):status='MIXED_PATHOLOGY_REPRESENTATION'
 else:status='FAIL_PATHOLOGY_REPRESENTATION'
 dump(A/'logs/GATE_RESULT.json',{'status':status,'panel_rules':panelrows,'negative_controls':negative,'C4B_allowed':status=='PASS_PATHOLOGY_REPRESENTATION','fits':len(list((A/'logs/runs').glob('*.json'))),'evaluations':len(p)})
 pd.DataFrame(panelrows).to_csv(A/'C4A_GAIN_SUMMARY.csv',index=False)
 runs=[json.loads(f.read_text()) for f in (A/'logs/runs').glob('*.json')];pd.DataFrame(runs).to_csv(A/'C4A_RESOURCE_PROFILE.csv',index=False)
 extraction=[json.loads(f.read_text()) for f in (A/'logs/extraction').glob('*.json')];rows=[]
 for e in CFG['encoders']:
  for r in extraction:
   sample=r['sample_id'];path=OUT/'embeddings'/e/f'{sample}.npz';rows.append({'encoder':e,'sample_id':sample,'dimension':768 if e=='pathology_encoder_1' else 512,'cache_path':str(path.relative_to(ROOT)),'cache_sha256':sha(path),'bytes':path.stat().st_size,'new_extraction_seconds':r['seconds'] if e=='pathology_encoder_1' else 0,'device':'CPU','order_verified':True,'estimated_FOV_um':r['estimated_FOV_um'],'source_h5_sha256':r['source_h5_sha256']})
 pd.DataFrame(rows).to_csv(A/'C4A_ENCODER_MANIFEST.csv',index=False)
 events=[json.loads(f.read_text()) for f in (A/'logs/registry').glob('*.json')]
 for r in events:
  for k,v in list(r.items()):
   if isinstance(v,(dict,list)):r[k]=json.dumps(v,sort_keys=True)
 pd.DataFrame(events).to_csv(A/'C4A_REGISTRY_EXTENSIONS.csv',index=False);print(status,flush=True)
if __name__=='__main__':analyze()
