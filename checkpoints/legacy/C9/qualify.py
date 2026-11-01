"""C9A metadata-only cohort and gene qualification. Never parses matrix values."""
from contract import *
import gzip
from PIL import Image

def main():
 rows=[];sets=[];files={};intersection=[]
 for s in SAMPLES:
  d=RAW/'filtered_count_matrices'/(s+'_filtered_count_matrix');sp=RAW/'spatial'/(s+'_spatial')
  files[s]={k:str((d/n).relative_to(ROOT)) for k,n in [('matrix','matrix.mtx.gz'),('features','features.tsv.gz'),('barcodes','barcodes.tsv.gz')]}
  with (d/'features.tsv.gz').open() as f:features=[l.strip() for l in f]
  with (d/'barcodes.tsv.gz').open() as f:bc=np.array([l.strip() for l in f])
  assert len(set(bc))==len(bc);sets.append(set(features))
  # Barcodes are the published filtered matrix locations; annotation and nCount columns unused.
  meta=pd.read_csv(RAW/'metadata'/(s+'_metadata.csv'),usecols=['patientid','subtype']);assert set(meta.patientid.astype(str))=={s}
  pos=pd.read_csv(sp/'tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode');pos=pos.loc[bc];excluded=int((pos.in_tissue!=1).sum());pos=pos[pos.in_tissue==1];bc=pos.index.to_numpy(str)
  sf=json.loads((sp/'scalefactors_json.json').read_text());xy=pos[['pixel_col','pixel_row']].to_numpy(float);mpp=55/sf['spot_diameter_fullres'];xyum=xy*mpp
  # Source lattice: 100um horizontal center pitch, vertical pitch sqrt(3)/2*100um.
  lattice=np.column_stack([pos.array_col.to_numpy()*50,pos.array_row.to_numpy()*np.sqrt(3)*50]);design=np.column_stack([lattice,np.ones(len(pos))]);pixel_coef=np.linalg.lstsq(design,xy,rcond=None)[0];pixel_sv=np.linalg.svd(pixel_coef[:2],compute_uv=False);mpp=1/float(np.mean(pixel_sv));xyum=xy*mpp;coef=np.linalg.lstsq(design,xyum,rcond=None)[0];res=np.linalg.norm(design@coef-xyum,axis=1);sing=np.linalg.svd(coef[:2],compute_uv=False)
  assert res.max()<2 and max(abs(sing-1))<.03,(s,res.max(),sing)
  hires=Image.open(sp/'tissue_hires_image.png');px=xy*sf['tissue_hires_scalef'];assert (px>=0).all() and (px[:,0]<hires.width).all() and (px[:,1]<hires.height).all()
  np.savez_compressed(CACHE/(s+'_metadata.npz'),barcodes=bc,xy_px=xy,xy_um=xyum,mpp=mpp,hires_scale=sf['tissue_hires_scalef'])
  rows.append({'dataset_id':'zenodo.4739739','study':'Wu et al. 2021, Nature Genetics, doi:10.1038/s41588-021-00911-1','patient_id':s,'sample_id':s,'disease':'breast cancer '+','.join(sorted(set(meta.subtype.astype(str)))),'platform':'Visium fresh frozen','num_locations':len(bc),'off_tissue_excluded':excluded,'num_assay_genes':len(features),'overlap_S2_50':len(set(parent.ORIGINAL['S2-50'])&set(features)),'overlap_S2_280':len(set(parent.ORIGINAL['S2-280'])&set(features)),'image_resolution':'pending PDF inspection','spatial_coordinates_available':True,'duplicate_check_pass':True,'cohort_role':'calibration' if s in CAL else 'evaluation','qualification_status':'PENDING_IMAGE_QUALIFICATION','exclusion_reason':'','mpp_fullres':mpp,'mpp_from_spot_diameter':55/sf['spot_diameter_fullres'],'source_lattice_max_residual_um':float(res.max()),'source_lattice_singular_values':sing.tolist(),'metadata_sha256':sha(CACHE/(s+'_metadata.npz'))})
 common=set.intersection(*sets);gg=[g for g in parent.ORIGINAL['S2-50'] if g in common];assert len(gg)>=40
 for p in ['S2-50','S2-280']:
  for g in parent.ORIGINAL[p]:intersection.append({'panel':p,'gene':g,'included':g in common,'reason':'exact symbol in all six assays' if g in common else 'absent from at least one assay'})
 pd.DataFrame(intersection).to_csv(A/'C9A_GENE_INTERSECTION.csv',index=False)
 dump(A/'logs/METADATA_QUALIFICATION.json',{'patients':rows,'genes':gg,'count_files':files,'expression_matrix_values_accessed':False,'clinical_annotations_not_used_as_predictors':True,'prior_metadata_preview':'Two metadata rows containing published nCount/nFeature values were incidentally displayed during file-format inspection. No per-gene values, reconstruction metrics, or outcomes were accessed; split and settings were not selected from these summaries.','sources':['https://zenodo.org/records/4739739','https://pmc.ncbi.nlm.nih.gov/articles/PMC9044823/'],'freshness':'Six source patient IDs in an Australian/independent-laboratory Visium cohort; distinct from four development Xenium patients, Swedish HER2ST A-H, GSE311609 diagnostic pilot and 10x Prime diagnostic release. No C0-C8 run used this study or these patients; published study includes separate HER2ST analyses, but those imported samples are NOT in this six-patient local cohort.'})
 print([(r['sample_id'],r['num_locations'],r['overlap_S2_50']) for r in rows],flush=True)
if __name__=='__main__':main()
