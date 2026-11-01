"""Image-only local cross-correlation landmarks and held-out geometric audit."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C7';sys.path.insert(0,str(OUT/'deps'))
import cv2,numpy as np
from PIL import Image
cv2.setNumThreads(4);cv2.setRNGSeed(0);dest=OUT/'C7A/registration'
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--base',choices=['shape','ecc'],default='shape');args=parser.parse_args();suffix='refined' if args.base=='shape' else 'ecc_audit'
for sample in ['GSM9509172','GSM9509160','GSM9509170','GSM9509165']:
 p=dest/(sample+'_'+args.base+'.json')
 if not p.exists():continue
 r=json.loads(p.read_text());f=np.load(OUT/'cache'/(sample+'_registration_arrays.npz'));d=f['dapi'];h=f['he'];mat=np.array(r['matrix_2um']);w=cv2.warpAffine(d,mat,h.shape[::-1]);hm=cv2.resize(f['mask_he'],h.shape[::-1],interpolation=cv2.INTER_NEAREST).astype(bool)
 hs=cv2.GaussianBlur(h.astype('float32'),(0,0),2);ws=cv2.GaussianBlur(w.astype('float32'),(0,0),2);pts=[]
 half=40;search=40
 for y in range(half+search,h.shape[0]-half-search,100):
  for x in range(half+search,h.shape[1]-half-search,100):
   if not hm[y,x]:continue
   patch=ws[y-half:y+half,x-half:x+half];target=hs[y-half-search:y+half+search,x-half-search:x+half+search]
   if patch.std()<5 or np.mean(w[y-half:y+half,x-half:x+half]>15)<.25:continue
   cc=cv2.matchTemplate(target,patch,cv2.TM_CCOEFF_NORMED);_,val,_,pos=cv2.minMaxLoc(cc);shift=np.array(pos)-search
   if val>.4 and np.max(abs(shift))<search-2:pts.append([x,y,x+shift[0],y+shift[1],val])
 pts=np.array(pts);rng=np.random.default_rng(0);order=rng.permutation(len(pts));test=order[::3];train=np.setdiff1d(order,test);cor,ins=cv2.estimateAffine2D(pts[train,:2],pts[train,2:4],method=cv2.RANSAC,ransacReprojThreshold=8,maxIters=20000,confidence=.999,refineIters=20)
 pred=pts[:,:2]@cor[:,:2].T+cor[:,2];err=np.linalg.norm(pred-pts[:,2:4],axis=1)*2
 final=(np.vstack([cor,[0,0,1]])@np.vstack([mat,[0,0,1]]))[:2];m=final.copy();m[:,:2]/=r['he_mpp'];m[:,2]*=2/r['he_mpp'];w2=cv2.warpAffine(d,final,h.shape[::-1]);im=Image.fromarray(np.stack([h,w2,h//2],-1));im.thumbnail((2200,1600));im.save(dest/(sample+'_'+suffix+'_overlay.png'))
 audit={**r,'dapi_um_to_he_pixel':m.tolist(),'matrix_2um':final.tolist(),'n_landmarks':len(pts),'train_inliers':int(ins.sum()),'heldout_landmarks':len(test),'heldout_median_error_um':float(np.median(err[test])),'heldout_p90_error_um':float(np.quantile(err[test],.9)),'heldout_fraction_error_le25um':float(np.mean(err[test]<=25)),'median_NCC':float(np.median(pts[:,4])),'status':'PENDING_VISUAL_REVIEW'}
 (dest/(sample+'_'+suffix+'.json')).write_text(json.dumps(audit,indent=2)+'\n');np.savez_compressed(dest/(sample+'_'+suffix+'_NCC_landmarks.npz'),landmarks_2um=pts,train=train,test=test,error_um=err);print(sample,'landmarks',len(pts),'test',len(test),'med',np.median(err[test]),'p90',np.quantile(err[test],.9),'fraction25',np.mean(err[test]<=25),'NCC',np.median(pts[:,4]),flush=True)
