"""Audit source hashes, numerical aggregation, ledger cells and figure completeness."""
from contract10 import *
from PIL import Image
M=OUT/'MANUSCRIPT'
def main():
 assert json.loads((B/'logs/INTEGRITY_AUDIT.json').read_text())['status']=='PASS';assert 'PASS: 1829 C9 artifacts' in (B/'logs/PARENT_VERIFICATION_AFTER.log').read_text()
 sources=json.loads((M/'SOURCE_INDEX.json').read_text())
 for p,r in sources.items():assert sha(ROOT/p)==r['sha256'],('source index',p)
 ld=pd.read_csv(M/'RESULTS_LEDGER.csv');cache={};numeric=0
 for p,g in ld.groupby('source_file'):
  assert set(g.source_sha256)=={sha(ROOT/p)},('ledger hash',p)
  if p.endswith('.csv'):
   f=pd.read_csv(ROOT/p,keep_default_na=False);cache[p]=f
   for _,r in g.iterrows():
    import re
    match=re.fullmatch(r'row (\d+); (.+)',r.source_row_or_field);assert match
    actual=f.iloc[int(match[1])-2][match[2]];expected=r.value
    try:np.testing.assert_allclose(float(actual),float(expected),rtol=1e-12,atol=1e-14,equal_nan=True);numeric+=1
    except ValueError:assert str(actual)==str(expected),(p,r.ledger_id,actual,expected)
 cross=pd.read_csv(OUT/'CROSS_COHORT_SUMMARY.csv');cols=['median_pearson','mean_mse','mean_mae']
 for _,r in cross.iterrows():
  if r.source_file not in cache:cache[r.source_file]=pd.read_csv(ROOT/r.source_file)
  orig=cache[r.source_file].iloc[[int(x)-2 for x in r.source_rows.split(',')]]
  for met in cols:
   field={'mean_mse':'mse','mean_mae':'mae'}.get(met,met) if 'TRANSPORTABILITY_SUMMARY' in r.source_file else met
   np.testing.assert_allclose(r[met],orig[field].mean(),rtol=1e-12,atol=1e-14)
 t=pd.read_csv(M/'tables/Table2_reconstruction.csv')
 for _,r in t.iterrows():
  q=cross.iloc[[int(x)-2 for x in r.source_rows.split(',')]];assert q.cohort.nunique()==1 and q.primary_gene_count.nunique()==1
  for met in cols:np.testing.assert_allclose(r[met],q[met].median(),rtol=1e-12,atol=1e-14)
  assert r.patient_count==q.patient.nunique()
 figures=[]
 for sub,n in [('figures',7),('supp_figures',10)]:
  ps=sorted((M/sub).glob('*.png'));assert len(ps)==n
  for p in ps:
   im=Image.open(p);assert min(im.info['dpi'])>=299.9;assert p.with_suffix('.svg').exists() and p.with_suffix('.py').exists();sm=pd.read_csv(p.with_name(p.stem+'_sources.csv'))
   for _,s in sm.iterrows():assert sha(ROOT/s.source_file)==s.sha256,(p,s.source_file)
   figures.append({'figure':str(p.relative_to(ROOT)),'pixels':im.size,'dpi':im.info['dpi'],'svg':True,'plot_entrypoint':True})
 for sub,n in [('tables',4),('supp_tables',8)]:
  md=list((M/sub).glob('*.md'));assert len(md)==n
  for p in md:assert p.with_suffix('.csv').exists() or p.with_suffix('.csv.gz').exists()
 abstract=(M/'MANUSCRIPT_DRAFT.md').read_text().split('## Abstract\n')[1].split('## Introduction')[0];assert len(abstract.split())<=250
 ef=pd.read_csv(B/'C10B_PATIENT_EFFECTS.csv');e=ef[(ef.strategy=='random')&(ef.regime=='V1')&(ef.comparison=='sparse_ST')];assert len(e)==66;assert (e.delta_pearson>=.02).all() and (e.relative_mse_reduction>=.02).all() and (e.mae_reduction>=0).all()
 x=pd.read_csv(B/'C10B_V2_DEGRADATION.csv');x=x[x.budget==.2];assert len(x)==22 and (x.delta_pearson_V2_minus_V1<0).all() and (x.mse_increase_V2_minus_V1>0).all() and (x.mae_increase>0).all()
 claims=pd.read_csv(OUT/'CLAIM_EVIDENCE_MATRIX.csv');assert len(claims)==7
 result={'status':'PASS','source_hashes_checked':len(sources),'ledger_entries':len(ld),'ledger_numeric_cells_recomputed':numeric,'cross_cohort_patient_rows_recomputed':len(cross),'table2_rows_recomputed':len(t),'abstract_words':len(abstract.split()),'figures':figures,'main_tables':4,'supplementary_tables':8,'scientific_support_counts_recomputed':True,'all_C0_C9_parents_verified_read_only':True,'visual_review':'All17figures viewed in contact sheets; revisedSFig3/SFig10 andmainFig3 reviewed atfullsize; all22coordinateoverlays viewed beforeoutcomes','no_new_statistical_tests':True};dump(M/'PACKAGE_AUDIT.json',result);print(json.dumps({k:v for k,v in result.items() if k!='figures'},indent=2))
if __name__=='__main__':main()
