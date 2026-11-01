"""Frozen ResNet18 and C4 heads: external H&E only, no expression access."""
from support import *
from PIL import Image
from torchvision.models import resnet18,ResNet18_Weights
import h5py
from threadpoolctl import threadpool_limits

def main():
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True);net=resnet18(weights=None);cp=ROOT/'data/weights/hub/checkpoints/resnet18-f37072fd.pth';net.load_state_dict(torch.load(cp,map_location='cpu',weights_only=True));net.fc=torch.nn.Identity();net.eval();transform=ResNet18_Weights.IMAGENET1K_V1.transforms()
 for s in SAMPLES:
  t=time.perf_counter();p=CACHE/(s+'_features.npz');patch=CACHE/(s+'_patches.h5')
  if not p.exists():
   vec=[];patch_hashes=[]
   with core.PeakRSS() as mem,h5py.File(patch) as f:
    bc=f['barcodes'][:].astype(str);xy=f['xy_um'][:]
    for k in range(0,len(bc),64):
     imgs=f['img'][k:k+64];patch_hashes.extend(hashlib.sha256(a.tobytes()).hexdigest() for a in imgs)
     with torch.inference_mode():vec.append(net(torch.stack([transform(Image.fromarray(a)) for a in imgs])).numpy())
    x=np.concatenate(vec);assert x.shape==(len(bc),512) and np.isfinite(x).all();np.savez_compressed(p,embeddings=x,barcodes=bc,xy_um=xy,patch_sha256=np.array(patch_hashes))
   dump(A/'logs'/f'{s}_EXTRACTION.json',{'sample_id':s,'seconds':time.perf_counter()-t,'peak_rss_mb':mem.peak/1024**2,'feature_sha256':sha(p),'patch_sha256':sha(patch),'encoder_checkpoint_sha256':sha(cp),'external_expression_accessed':False,'source_commit':commit()})
  with np.load(p) as f:x=f['embeddings'].copy();bc=f['barcodes'].copy()
  for panel,genes in panels().items():
   dst=CACHE/(s+'_'+panel+'_image.npz')
   if dst.exists():continue
   rows=[json.loads(p.read_text()) for p in (ROOT/'C4/C4A/logs/runs').glob('*.json')];rows=[r for r in rows if r['encoder']==ENC and r['panel']==panel and r['model']=='D1' and r['phase']=='primary'];assert len(rows)==3;ys=[];refs=[]
   for r in sorted(rows,key=lambda r:r['seed']):
    checkpoint=ROOT/r['model_path'];assert sha(checkpoint)==r['model_sha256'];model=core.mlp_model(r['input_dim'],r['output_dim'],r['params']);model.load_state_dict(torch.load(checkpoint,map_location='cpu',weights_only=True));model.eval()
    with np.load(ROOT/r['preprocessing']/'feature.npz') as f:xx=((x-f['mean'])/f['scale']).astype(np.float32)
    with np.load(ROOT/r['preprocessing']/'target.npz') as f:y=core.predict(model,xx,'B2').astype(float)*f['scale']+f['mean']
    ys.append(y);refs.append({'checkpoint':r['model_path'],'sha256':r['model_sha256'],'preprocessing':r['preprocessing'],'feature_scaler_sha256':sha(ROOT/r['preprocessing']/'feature.npz'),'target_scaler_sha256':sha(ROOT/r['preprocessing']/'target.npz')})
   out=np.mean(ys,axis=0).astype(np.float32)[:,[ORIGINAL[panel].index(g) for g in genes]];np.savez_compressed(dst,prediction=out,genes=np.array(genes),barcodes=bc);dump(dst.with_suffix('.json'),{'sha256':sha(dst),'parents':refs,'external_labels_used':False})
  print('EXTRACTED',s,len(x),round(time.perf_counter()-t,2),flush=True)
if __name__=='__main__':
 with threadpool_limits(4):main()
