"""Version copies of frozen C0 masks without expression access."""
from common import *
import shutil
if __name__=='__main__':
 rows=[]
 for r in json.loads((ROOT/'reports/C0/C0_MASK_INDEX.json').read_text()):
  source=ROOT/r['path'];assert sha(source)==r['sha256'];dest=OUT/'masks'/source.name;shutil.copyfile(source,dest)
  with np.load(dest) as f:obs=f['observed'];bar=f['barcodes'];dist=f['distance_to_anchor_um'];assert obs.dtype==bool and obs.any() and (~obs).any()
  with np.load(ROOT/'C4/embeddings/resnet18_reference'/f'{r["sample_id"]}.npz') as f:assert np.array_equal(bar,f['barcodes'])
  row={'mask_id':r['mask_id'],'sample_id':r['sample_id'],'patient_id':SPLITS['patients'][r['sample_id']],'mask_type':{'random':'M1','block':'M2','capture_window':'M3'}[r['regime']],'regime':r['regime'],'nominal_coverage':r['requested_spot_coverage'],'realized_coverage':float(obs.mean()),'num_observed':int(obs.sum()),'num_hidden':int((~obs).sum()),'geometry':r['geometry'],'seed':r['seed'],'replicate':r['replicate'],'checksum':sha(dest),'path':str(dest.relative_to(ROOT)),'C0_source_path':r['path'],'center_um':r.get('center_um'),'radius_um':r.get('radius_um')};rows.append(row)
 dump(OUT/'configs/masks.json',rows);pd.DataFrame(rows).to_csv(OUT/'C5_MASK_INDEX.csv',index=False);print('Verified/copied',len(rows),'frozen masks')
