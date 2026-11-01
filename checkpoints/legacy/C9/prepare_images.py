"""Frozen H&E feature/head inference, without any expression matrix access."""
from contract import *
from PIL import Image
import torch
from torchvision.models import resnet18,ResNet18_Weights
from threadpoolctl import threadpool_limits

def main():
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True);net=resnet18(weights=None);cp=ROOT/'data/weights/hub/checkpoints/resnet18-f37072fd.pth';net.load_state_dict(torch.load(cp,map_location='cpu',weights_only=True));net.fc=torch.nn.Identity();net.eval();transform=ResNet18_Weights.IMAGENET1K_V1.transforms();q=json.loads((A/'logs/METADATA_QUALIFICATION.json').read_text());gg=q['genes'];info={}
 for s in SAMPLES:
  start=time.perf_counter();sp=RAW/'spatial'/(s+'_spatial');path=sp/'tissue_hires_image.png';im=Image.open(path).convert('RGB')
  with np.load(CACHE/(s+'_metadata.npz')) as f:bc=f['barcodes'];xy=f['xy_px']*f['hires_scale'];xyum=f['xy_um'];mpp=float(f['mpp']/f['hires_scale'])
  width=parent.count(112/mpp,1);vec=[];patch_hashes=[]
  with core.PeakRSS() as mem:
   for k in range(0,len(bc),64):
    patches=[]
    for x,y in xy[k:k+64]:
     left=int(np.floor(x-width/2+.5));top=int(np.floor(y-width/2+.5));a=im.crop((left,top,left+width,top+width)).resize((224,224),Image.Resampling.LANCZOS);patch_hashes.append(hashlib.sha256(np.array(a).tobytes()).hexdigest());patches.append(transform(a))
    with torch.inference_mode():vec.append(net(torch.stack(patches)).numpy())
   x=np.concatenate(vec);fp=CACHE/(s+'_features.npz');np.savez_compressed(fp,embeddings=x,barcodes=bc,xy_um=xyum,patch_sha256=patch_hashes)
   rows=[json.loads(p.read_text()) for p in (ROOT/'C4/C4A/logs/runs').glob('*.json')];rows=[r for r in rows if r['encoder']==ENC and r['panel']=='S2-50' and r['model']=='D1' and r['phase']=='primary'];assert len(rows)==3;ys=[];refs=[]
   for r in sorted(rows,key=lambda r:r['seed']):
    checkpoint=ROOT/r['model_path'];assert sha(checkpoint)==r['model_sha256'];model=core.mlp_model(r['input_dim'],r['output_dim'],r['params']);model.load_state_dict(torch.load(checkpoint,map_location='cpu',weights_only=True));model.eval()
    with np.load(ROOT/r['preprocessing']/'feature.npz') as f:xx=((x-f['mean'])/f['scale']).astype(np.float32)
    with np.load(ROOT/r['preprocessing']/'target.npz') as f:y=core.predict(model,xx,'B2').astype(float)*f['scale']+f['mean']
    ys.append(y);refs.append({'checkpoint':r['model_path'],'sha256':sha(checkpoint),'feature_scaler':str(Path(r['preprocessing'])/'feature.npz'),'target_scaler':str(Path(r['preprocessing'])/'target.npz')})
   pred=np.mean(ys,axis=0).astype(np.float32)[:,[parent.ORIGINAL['S2-50'].index(g) for g in gg]];np.savez_compressed(CACHE/(s+'_image.npz'),prediction=pred,genes=gg,barcodes=bc)
  rec={'sample_id':s,'feature_sha256':sha(fp),'image_prediction_sha256':sha(CACHE/(s+'_image.npz')),'encoder_sha256':sha(cp),'source_image_sha256':sha(path),'native_patch_pixels':width,'mpp':mpp,'nominal_FOV_um':112,'actual_FOV_um':width*mpp,'no_expression_values_accessed':True,'seconds':time.perf_counter()-start,'peak_rss_mb':mem.peak/1024**2,'parents':refs};dump(A/'logs'/(s+'_EXTRACTION.json'),rec);info[s]={'resolution':f'{im.width}x{im.height}, {mpp:.4f} um/pixel; {width} native pixels per112um patch','mpp':mpp,'native_patch_pixels':width};print('IMAGE',s,len(bc),round(rec['seconds'],1),flush=True)
 dump(A/'logs/IMAGE_QUALIFICATION.json',{'status':'LIMITED_FRESH_EXTERNAL_DATASET','input_contract':'Supplied unannotated Space Ranger tissue_hires_image.png, source-scaled coordinates, 112um square using lattice-calibrated mpp, nearest-integer width/top-left, LANCZOS to224 RGB; exact C0 ResNet18 IMAGENET1K_V1 transform and frozen weights.','patients':info,'limitation':'Only15-25 native pixels per112um patch in supplied common images. Image-resolution domain shift weakens M0 and limits H&E comparisons; no superresolution or encoder change. PDF inspection finds first4patients similarly2000px whole-image views and larger crops for2patients only, without common fullres coordinate transforms. No PDF annotations used as model input. Molecular comparisons remain on original registered spots.','PDF_inspected':True,'image_outcome_selection':False})
if __name__=='__main__':
 with threadpool_limits(4):main()
