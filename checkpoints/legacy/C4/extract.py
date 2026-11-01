"""Frozen public Phikon CLS features, official preprocessing, exact C0 patch identities."""
from common import *
sys.path.insert(0,str(OUT/'deps'))
import argparse,time,h5py,torch
from transformers import ViTImageProcessor,ViTModel

def setup():
 torch.set_num_threads(8);torch.use_deterministic_algorithms(True)
 folder=OUT/'weights/phikon';processor=ViTImageProcessor.from_pretrained(folder,local_files_only=True)
 model=ViTModel.from_pretrained(folder,local_files_only=True,add_pooling_layer=False,attn_implementation='sdpa').eval()
 for p in model.parameters():p.requires_grad_(False)
 return processor,model

def audit(sample):
 rows={r['sample_id']:r for r in json.loads((ROOT/'reports/C0/C0_DATA_MANIFEST.json').read_text())['slides']};r=rows[sample]
 path=ROOT/'data/raw/IDC/patches'/f'{sample}.h5'
 with h5py.File(path) as f,h5py.File(ROOT/'data/raw/IDC/adata'/f'{sample}.h5ad') as ad:
  barcodes=f['barcode'][:].reshape(-1).astype(str);lookup={b:i for i,b in enumerate(ad['obs/_index'][:].astype(str))};xy=ad['obsm/spatial'][:][[lookup[b] for b in barcodes]];error=float(np.max(abs(f['coords'][:][:,::-1]-xy)));native=int(f['img'].attrs['patch_size'])
  assert f['img'].shape==(len(barcodes),224,224,3) and f['img'].dtype==np.uint8 and error<=1.01
 with np.load(ROOT/'data/cache'/f'{sample}_features.npz') as f:assert np.array_equal(barcodes,f['barcodes'])
 with np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as f:xy_um=f['xy_um'];assert np.array_equal(barcodes,f['barcodes'])
 return barcodes,xy_um,{'sample_id':sample,'patient':r['patient'],'n_locations':len(barcodes),'native_patch_pixels':native,'C0_estimated_native_mpp':r['microns_per_pixel'],'estimated_FOV_um':native*r['microns_per_pixel'],'estimated_resized_mpp':native*r['microns_per_pixel']/224,'alignment_max_error_px':error,'source_h5_sha256':sha(path),'RGB_shape':[224,224,3],'C4_patch_geometry_changed':False}

def extract(sample,processor,model):
 folder=OUT/'embeddings/pathology_encoder_1';folder.mkdir(exist_ok=True);path=folder/f'{sample}.npz';meta=A/'logs/extraction'/f'{sample}.json'
 if meta.exists():assert sha(path)==json.loads(meta.read_text())['cache_sha256'];return
 start=time.perf_counter();barcodes,xy,record=audit(sample);vectors=[];hashes=[]
 with core.PeakRSS() as memory,h5py.File(ROOT/'data/raw/IDC/patches'/f'{sample}.h5') as f:
  for begin in range(0,len(barcodes),16):
   images=f['img'][begin:begin+16];hashes.extend(hashlib.sha256(image.tobytes()).hexdigest() for image in images)
   x=processor(list(images),return_tensors='pt').pixel_values
   with torch.inference_mode():vectors.append(model(x).last_hidden_state[:,0,:].numpy())
   if begin%512==0:print(sample,begin,'/',len(barcodes),'elapsed',round(time.perf_counter()-start),flush=True)
  embeddings=np.concatenate(vectors);assert embeddings.shape==(len(barcodes),768) and np.isfinite(embeddings).all()
  enc=json.loads((A/'configs/encoder_files.json').read_text());phash=sha(A/'configs/extraction.json')
  np.savez_compressed(path,embeddings=embeddings,barcodes=barcodes,xy_um=xy,source_patch_sha256=np.array(hashes),source_patch_row=np.arange(len(barcodes)),sample_id=np.array(sample),patient_id=np.array(record['patient']),encoder_preprocessing_sha256=np.array(phash))
  # New reference wrapper contains the exact old feature values; C0 is never overwritten or extracted again.
  e0=OUT/'embeddings/resnet18_reference';e0.mkdir(exist_ok=True)
  with np.load(ROOT/'data/cache'/f'{sample}_features.npz') as old:x0=old['embeddings'].copy()
  np.savez_compressed(e0/f'{sample}.npz',embeddings=x0,barcodes=barcodes,xy_um=xy,source_patch_sha256=np.array(hashes),source_patch_row=np.arange(len(barcodes)),sample_id=np.array(sample),patient_id=np.array(record['patient']),encoder_preprocessing_sha256=np.array(sha(ROOT/'data/cache'/f'{sample}_features.json')))
 r={**record,'encoder':'phikon','encoder_revision':enc['revision'],'embedding_dimension':768,'seconds':time.perf_counter()-start,'peak_rss_mb':memory.peak/1024**2,'threads':8,'device':'CPU','precision':'float32','batch_size':16,'cache_path':str(path.relative_to(ROOT)),'cache_sha256':sha(path),'cache_bytes':path.stat().st_size,'reference_cache_sha256':sha(e0/f'{sample}.npz'),'preprocessing_sha256':phash,'source_commit':commit()}
 dump(meta,r);registry('extraction',f'C4A_extract_phikon_{sample}',r);print('DONE',sample,r['seconds'],flush=True)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['trainval','test']);args=parser.parse_args()
 if args.stage=='test':assert (A/'SELECTION_LOCK.json').exists(),'Test extraction follows selection lock'
 processor,model=setup()
 for sample in SPLITS['splits']['test'] if args.stage=='test' else SPLITS['splits']['train']+SPLITS['splits']['validation']:extract(sample,processor,model)
if __name__=='__main__':main()
