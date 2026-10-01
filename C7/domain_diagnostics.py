"""Post-evaluation descriptive domain-shift diagnostics; no calibration or fitting."""
from support import *
def main():
 rows=[]
 for sample in SAMPLES:
  for panel in panels():
   y=truth(sample,panel);image=scene(sample,panel)['image'];bias=image.mean(0)-y.mean(0)
   with np.load(CACHE/(sample+'_targets.npz')) as t:
    idx=[list(t['genes']).index(g) for g in panels()[panel]];raw=t['raw_counts'][:,idx]
   rows.append({'sample_id':sample,'panel':panel,'zero_fraction':float((raw==0).mean()),'target_mean':float(y.mean()),'image_prediction_mean':float(image.mean()),'mean_gene_target_variance':float(y.var(0).mean()),'full_tissue_E0_mse':float(((y-image)**2).mean()),'mean_gene_squared_mean_bias':float((bias**2).mean())})
 pd.DataFrame(rows).to_csv(B/'C7B_DOMAIN_SHIFT_DIAGNOSTICS.csv',index=False)
if __name__=='__main__':main()
