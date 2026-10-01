"""Complete provenance-rich tables, distance/support summaries and audit inventory."""
from synthesis import *
import platform

def main():
 # Preserve existing source index across separate package stages.
 SOURCES.update(json.loads((M/'SOURCE_INDEX.json').read_text()))
 c0=read('reports/C0/C0_DATASET_STATISTICS.csv');h=read('C7/C7A_EXTERNAL_DATASET_MANIFEST.csv');h=h[h.dataset_id=='HER2ST'];w=read('C9/C9A/C9A_EXTERNAL_DATASET_MANIFEST.csv');b=read('C10/C10A_DATASET_MANIFEST.csv')
 records=[{'cohort':'Development','study':'C0 selected breast datasets; source URLs in C0 statistics','disease':'Breast cancer','platform':'Xenium, 280 shared targeted genes','patient_count':4,'evaluation_patient_count':1,'location_count':int(c0.n_locations.sum()),'evaluation_locations':int(c0[c0.split=='test'].n_locations.sum()),'H&E resolution':'0.2125–0.3639 um/pixel; HEST estimated physical scale','primary_gene_overlap':'50/50','role':'2 train,1 validation,1 primary test;4 conditional LOPO sensitivity','source_file':'reports/C0/C0_DATASET_STATISTICS.csv','source_rows':'2–5'},
 {'cohort':'HER2ST','study':'Andersson et al.; doi:10.1038/s41467-021-26271-2','disease':'HER2-positive breast cancer','platform':'Spatial Transcriptomics,100umspots/200umpitch','patient_count':8,'evaluation_patient_count':8,'location_count':int(h.num_locations.sum()),'evaluation_locations':int(h.num_locations.sum()),'H&E resolution':f'{h.image_resolution.min():.4f}–{h.image_resolution.max():.4f} um/pixel','primary_gene_overlap':'40/50 common','role':'C7 zero-shot; C8 transportability development. Not fresh after C7 inspection.','source_file':'C7/C7A_EXTERNAL_DATASET_MANIFEST.csv','source_rows':','.join(map(str,h.source_row))},
 {'cohort':'Wu fresh validation','study':w.study.iloc[0],'disease':'Breast cancer; local ER-positive and TNBC','platform':'Visium fresh frozen','patient_count':len(w),'evaluation_patient_count':int((w.cohort_role=='evaluation').sum()),'location_count':int(w.num_locations.sum()),'evaluation_locations':int(w[w.cohort_role=='evaluation'].num_locations.sum()),'H&E resolution':'15/25 native pixels per112um patch;2000px previews','primary_gene_overlap':'50/50','role':'3 calibration;3 heldout evaluation','source_file':'C9/C9A/C9A_EXTERNAL_DATASET_MANIFEST.csv','source_rows':','.join(map(str,w.source_row))},
 {'cohort':'Bassiouni second replication','study':b.study.iloc[0],'disease':'Primary TNBC','platform':'Visium fresh frozen','patient_count':len(b),'evaluation_patient_count':len(b),'location_count':int(b.num_locations.sum()),'evaluation_locations':int(b.num_locations.sum()),'H&E resolution':'24–29 native pixels per112um patch;2000px previews','primary_gene_overlap':'50/50','role':'22 evaluation;first section perpatient;no C10 calibration','source_file':'C10/C10A_DATASET_MANIFEST.csv','source_rows':','.join(map(str,b.source_row))}]
 save(pd.DataFrame(records),'Table1_cohorts','tables',True)
 # Complete exact intersection table.
 panels=pd.DataFrame([{'checkpoint':'C0','panel':panel,'gene':g,'included':True} for panel,gg in parent.ORIGINAL.items() for g in gg]);parts=[panels];d=read('C10/C10A/GENE_INTERSECTION.csv');d['checkpoint']='C10';parts.append(d);SOURCES.pop('C10/MANUSCRIPT/supp_tables/SuppTable1_gene_panels.csv',None)
 for cp,path in [('C7','C7/C7A/C7A_GENE_INTERSECTION.csv'),('C9','C9/C9A/C9A_GENE_INTERSECTION.csv')]:
  d=read(path);d['checkpoint']=cp
  if 'included_in_evaluation' in d:d['included']=d.included_in_evaluation
  parts.append(d)
 save(pd.concat(parts,ignore_index=True),'SuppTable1_gene_panels','supp_tables',True)
 # Distance summaries use equal geometry means, then descriptive within-patient occupied-bin trends.
 dd=read('C10/C10B/C10B_DISTANCE_TO_ANCHOR.csv');q=dd[dd.control.isna()&(dd.strategy=='random')&(dd.budget==.2)&(dd.method=='M1')];keys=['patient_id','regime','lower_um'];per=q.groupby(keys+['geometry'])[['mse','mae','n_hidden','median_pearson']].mean().groupby(keys).mean().reset_index();per.to_csv(B/'C10B_DISTANCE_PATIENT_CURVES.csv',index=False);trend=[]
 for key,g in per.groupby(['patient_id','regime']):
  g=g.sort_values('lower_um');dx=np.diff(g.mse);trend.append({'patient_id':key[0],'regime':key[1],'occupied_bins':len(g),'strictly_increasing_mse':bool(len(dx)>0 and (dx>0).all()),'farthest_minus_nearest_mse':float(g.mse.iloc[-1]-g.mse.iloc[0]),'source_file':'C10/C10B/C10B_DISTANCE_PATIENT_CURVES.csv','source_rows':','.join(map(str,g.index+2))})
 trends=save(pd.DataFrame(trend),'C10_distance_trends');trends.to_csv(B/'C10B_DISTANCE_TRENDS.csv',index=False)
 effects=read('C10/C10B/C10B_PATIENT_EFFECTS.csv');primary=effects[(effects.regime=='V1')&(effects.strategy=='random')&(effects.budget==.2)];im=primary[primary.comparison=='image'];degr=read('C10/C10B/C10B_V2_DEGRADATION.csv');degr=degr[degr.budget==.2]
 summaries=[]
 def stat(name,value,unit,frame,field):summaries.append({'statistic':name,'value':value,'unit':unit,'source_file':'|'.join(frame.source_file.unique()),'source_rows':','.join(map(str,frame.source_row)),'source_field':field})
 stat('fusion_patients_Pearson_improved',int((im.delta_pearson>0).sum()),'patients',im,'count(delta_pearson>0)');stat('fusion_patients_MSE_improved',int((im.mse_reduction>0).sum()),'patients',im,'count(mse_reduction>0)');stat('fusion_patients_MAE_improved',int((im.mae_reduction>0).sum()),'patients',im,'count(mae_reduction>0)');stat('fusion_median_delta_Pearson',im.delta_pearson.median(),'Pearson difference',im,'median(delta_pearson)');stat('fusion_median_relative_MSE_increase',-im.relative_mse_reduction.median(),'fraction',im,'-median(relative_mse_reduction)');stat('V2_all_metrics_worse_patients',int(((degr.delta_pearson_V2_minus_V1<0)&(degr.mse_increase_V2_minus_V1>0)&(degr.mae_increase>0)).sum()),'patients',degr,'joint count');stat('V2_median_relative_MSE_increase',degr.relative_mse_increase.median(),'fraction',degr,'median(relative_mse_increase)')
 for reg in ['V1','V2']:
  tr=read('C10/C10B/C10B_DISTANCE_TRENDS.csv');tr=tr[tr.regime==reg];stat(reg+'_distance_strictly_increasing_patients',int(tr.strictly_increasing_mse.sum()),'patients',tr,'count(strictly_increasing_mse)');stat(reg+'_distance_farthest_worse_patients',int((tr.farthest_minus_nearest_mse>0).sum()),'patients',tr,'count(farthest_minus_nearest_mse>0)')
 save(pd.DataFrame(summaries),'C10_support_statistics')
 # Frozen secondary FPS results; no post-hoc strategy selection.
 pt=read('C10/C10B/C10B_PATIENT_METRICS.csv');fps=[]
 for key,g in pt[(pt.method=='M1')&(pt.budget==.2)].groupby(['patient_id','regime']):
  r=g.set_index('strategy');fps.append({'patient':key[0],'regime':key[1],'FPS_minus_random_pearson':r.loc['spatial_FPS','median_pearson']-r.loc['random','median_pearson'],'random_minus_FPS_MSE':r.loc['random','mean_mse']-r.loc['spatial_FPS','mean_mse'],'source_file':'C10/C10B/C10B_PATIENT_METRICS.csv','source_rows':','.join(map(str,g.source_row))})
 save(pd.DataFrame(fps),'C10_FPS_effects')
 audit=json.loads((B/'logs/INTEGRITY_AUDIT.json').read_text());old=json.loads((ROOT/'C9/C9_STATUS.json').read_text())['technical_audit'];save(pd.DataFrame([{'checkpoint':'C9',**old,'source_file':'C9/C9_STATUS.json','source_field':'technical_audit'},{'checkpoint':'C10',**audit,'source_file':'C10/C10B/logs/INTEGRITY_AUDIT.json','source_field':'root'}]),'SuppTable5_integrity','supp_tables',True)
 masks=[]
 for cp in ['C9','C10']:
  for p in (ROOT/cp/(cp+'B')/'masks').glob('*.json'):
   r=json.loads(p.read_text());masks.append({**r,'checkpoint':cp,'source_file':str(p.relative_to(ROOT)),'source_sha256':sha(p)})
 save(pd.DataFrame(masks),'SuppTable6_masks','supp_tables',True)
 timeline=[]
 statuses={0:'PASS_C0 (see C0 report)',1:'PASS_DEVELOPMENT (accepted)',2:'MIXED_CONTEXT_SIGNAL',3:'MIXED_SPATIAL_SPECIFICITY',4:'MIXED_SPATIAL_RECHECK',5:'PASS_SPARSE_ST_RECONSTRUCTION',6:'MIXED_SAMPLING_EFFECT',7:'MIXED_EXTERNAL_VALIDATION',8:'MIXED_CROSS_PLATFORM_CALIBRATION',9:'MIXED_FRESH_EXTERNAL_VALIDATION',10:'PASS_LIMITED_SECOND_REPLICATION'}
 for i in range(11):
  path=ROOT/f'C{i}/C{i}_STATUS.json';path=path if path.exists() else ROOT/('reports/C0/C0_REPORT.md' if i==0 else f'C{i}/C{i}_REPORT.md');timeline.append({'checkpoint':f'C{i}','descriptive_status':statuses[i],'source_file':str(path.relative_to(ROOT)),'source_sha256':sha(path),'exact_status_payload':path.read_text() if path.suffix=='.json' else 'See report; user accepted C1 as PASS_DEVELOPMENT'})
 save(pd.DataFrame(timeline),'SuppTable8_checkpoint_timeline','supp_tables',True)
 gate=json.loads((B/'logs/SCIENTIFIC_GATE.json').read_text());table4=[{'cohort':'Wu fresh validation','held_out_patients':3,'primary_gene_count':50,'sparse_ST_beats_HE':'3/3patients at all3budgets, r/MSE/MAE','coverage_response':'3/3patients','distance_effect':'V1 increasing occupied-bin error in3/3;V2 nonmonotonic','contiguous_extrapolation':'3/3worse at20%','HE_incremental_gain':'1/3Pearson,0/3MSE/MAE at20%','validation_status':'MIXED_FRESH_EXTERNAL_VALIDATION; core supported; affine effect not identifiable; image resolution limited','source_file':'C9/C9_REPORT.md;C9/C9B/C9B_PATIENT_EFFECTS.csv'},
 {'cohort':'Bassiouni second replication','held_out_patients':22,'primary_gene_count':50,'sparse_ST_beats_HE':'22/22at all3budgets meet joint meaningful criterion','coverage_response':str(gate['coverage_response_patients'])+'/22patients','distance_effect':str(trends.groupby('regime').strictly_increasing_mse.sum().to_dict())+' strictly increasing; see bins/counts','contiguous_extrapolation':'22/22worse r/MSE/MAE at20%','HE_incremental_gain':f'{(im.delta_pearson>0).sum()}/22Pearson;0/22MSE/MAE at20%','validation_status':gate['status_pending_technical_audit'],'source_file':'C10/C10B/logs/SCIENTIFIC_GATE.json;C10/C10B/C10B_PATIENT_EFFECTS.csv;C10/C10B/C10B_DISTANCE_TRENDS.csv'}];save(pd.DataFrame(table4),'Table4_external_replication','tables',True)
 dump(M/'SOURCE_INDEX.json',SOURCES)
 env={'python':sys.version,'platform':platform.platform(),'processor':platform.processor(),'numpy':np.__version__,'pandas':pd.__version__,'threads':4,'execution':'CPU; frozen image forward passes, interpolation and audits only; noC10training','lscpu':subprocess.check_output(['lscpu'],text=True)};dump(OUT/'C10B/logs/ENVIRONMENT.json',env)
 print('SUPPLEMENT COMPLETE',pd.DataFrame(summaries).to_string(index=False),flush=True)
if __name__=='__main__':main()
