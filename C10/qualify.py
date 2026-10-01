"""Qualify using only gene/barcode/image/coordinate metadata, never H5 count values."""
from contract10 import *
from PIL import Image

def main():
 samples=json.loads((A/'logs/SOURCE_SAMPLE_SELECTION.json').read_text())['selected'];inputs={};sets=[];rows=[]
 for r in samples:
  s=r['patient_id'];paths={}
  for key,suffix in [('counts','filtered_feature_bc_matrix.h5'),('image','tissue_hires_image.png'),('coordinates','tissue_positions_list.csv'),('scales','scalefactors_json.json')]:
   matches=list(RAW.glob(r['gsm']+'*'+suffix));assert len(matches)==1;paths[key]=str(matches[0].relative_to(ROOT))
  inputs[s]=paths
  with h5py.File(ROOT/paths['counts']) as f:gg=f['matrix/features/name'][:].astype(str);bc=f['matrix/barcodes'][:].astype(str);shape=f['matrix/shape'][:];assert shape[1]==len(bc) and shape[0]==len(gg)
  assert len(set(bc))==len(bc);sets.append(set(gg))
  pos=pd.read_csv(ROOT/paths['coordinates'],header=None,names=['barcode','in_tissue','row','col','y','x']).set_index('barcode').loc[bc];excluded=int((pos.in_tissue!=1).sum());pos=pos[pos.in_tissue==1];bc=pos.index.to_numpy(str);xy=pos[['x','y']].to_numpy(float)
  sf=json.loads((ROOT/paths['scales']).read_text());lattice=np.column_stack([pos.col.to_numpy()*50,pos.row.to_numpy()*np.sqrt(3)*50]);design=np.column_stack([lattice,np.ones(len(pos))]);coef=np.linalg.lstsq(design,xy,rcond=None)[0];sv=np.linalg.svd(coef[:2],compute_uv=False);mpp=1/np.mean(sv);res=np.linalg.norm(design@coef-xy,axis=1)*mpp;anisotropy=max(sv)/min(sv)
  assert res.max()<2 and anisotropy<1.03,(s,res.max(),anisotropy)
  im=Image.open(ROOT/paths['image']);px=xy*sf['tissue_hires_scalef'];assert (px>=0).all() and (px[:,0]<im.width).all() and (px[:,1]<im.height).all();assert len(bc)>=100
  np.savez_compressed(CACHE/(s+'_metadata.npz'),barcodes=bc,xy_px=xy,xy_um=xy*mpp,mpp=mpp,hires_scale=sf['tissue_hires_scalef'])
  mm=mpp/sf['tissue_hires_scalef'];width=parent.count(112/mm,1)
  rows.append({'dataset':'GSE210616','study':'Bassiouni et al. Cancer Research 2023; doi:10.1158/0008-5472.CAN-22-2682','patient_id':s,'sample_id':r['gsm'],'disease':'primary triple-negative breast cancer','platform':'Visium fresh frozen','num_locations':len(bc),'gene_overlap_S2_50':len(set(parent.ORIGINAL['S2-50'])&set(gg)),'gene_overlap_S2_280':len(set(parent.ORIGINAL['S2-280'])&set(gg)),'H&E resolution':f'{im.width}x{im.height}; {mm:.4f}um/pixel; {width}px/112um patch','coordinate availability':True,'prior-use check':'source identity and repository identifier/hash review','qualification status':'LIMITED_SECOND_REPLICATION_DATASET','reason':'Public H&E downsampled; frozen image baseline has resolution domain shift. All22patients used; no C10 target fitting.','off_tissue_excluded':excluded,'max_lattice_residual_um':float(res.max()),'anisotropy':float(anisotropy),'native_patch_pixels':width})
 common=set.intersection(*sets);gene=[g for g in parent.ORIGINAL['S2-50'] if g in common];assert len(gene)>=40
 for s in SAMPLES:
  with h5py.File(ROOT/inputs[s]['counts']) as f:gg=list(f['matrix/features/name'][:].astype(str));assert all(gg.count(g)==1 for g in gene)
 pd.DataFrame(rows).to_csv(OUT/'C10A_DATASET_MANIFEST.csv',index=False)
 pd.DataFrame([{'panel':panel,'gene':g,'included':g in common,'reason':'exact symbol in all assays' if g in common else 'absent'} for panel in ['S2-50','S2-280'] for g in parent.ORIGINAL[panel]]).to_csv(A/'GENE_INTERSECTION.csv',index=False)
 dump(A/'logs/QUALIFICATION.json',{'status':'LIMITED_SECOND_REPLICATION_DATASET','patients':rows,'inputs':inputs,'genes':gene,'expression_values_accessed':False,'freshness':'22 source-identified TNBC patients from Virginia Commonwealth University / USC study; not the original Xenium patients, Swedish HER2ST, Wu local Australian cohort, GSE311609 pilot or Prime diagnostic release','sources':['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616','https://doi.org/10.1158/0008-5472.CAN-22-2682'],'selection':'First numbered section for every published patient; no performance or expression-derived selection','source_ids_before_outcome':True})
 print([(r['patient_id'],r['num_locations'],r['gene_overlap_S2_50'],r['native_patch_pixels']) for r in rows],flush=True)
if __name__=='__main__':main()
