"""Describe rank sensitivity from floating-point ties without changing locked metrics."""
from contract import *
def main():
 rows=[]
 for path in sorted((B/'logs/runs').glob('*.json')):
  r=json.loads(path.read_text())
  with np.load(ROOT/r['prediction_path']) as f:
   delta=f['M2_spearman']-f['M1_spearman'];finite=delta[np.isfinite(delta)]
   rows.append({'run_id':r['run_id'],'patient_id':r['patient_id'],'regime':r['regime'],'budget':r['budget'],'control':r['control'],'prediction_max_abs_difference':float(np.max(abs(f['M2']-f['M1']))),'max_abs_gene_spearman_difference':float(np.max(abs(finite))) if len(finite) else 0.,'mean_gene_spearman_difference':float(np.mean(finite)) if len(finite) else 0.})
 pd.DataFrame(rows).to_csv(B/'C9B_NUMERICAL_RANK_SENSITIVITY.csv',index=False)
 dump(B/'logs/NUMERICAL_AUDIT.json',{'status':'DISCLOSED_NUMERICAL_LIMITATION','maximum_prediction_difference':max(r['prediction_max_abs_difference'] for r in rows),'maximum_gene_spearman_difference':max(r['max_abs_gene_spearman_difference'] for r in rows),'cause':'Raw Spearman rankdata breaks machine-precision ties differently after affine transform/inverse. Pearson and common-unit errors retain numerical equivalence. Raw locked metrics are preserved; rank differences are not harmonization evidence.','no_prediction_or_metric_replacement':True,'no_status_or_selection_change':True})
if __name__=='__main__':main()
