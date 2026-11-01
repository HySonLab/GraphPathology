"""Label-free H&E/DAPI affine registration; expression is never opened."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C7';sys.path.insert(0,str(OUT/'deps'))
import cv2,numpy as np,tifffile,openslide
from skimage.color import rgb2hed
from PIL import Image
cv2.setNumThreads(4);cv2.setRNGSeed(0)
DEST=OUT/'C7A/registration';DEST.mkdir(exist_ok=True)
def norm(a):
 lo,hi=np.percentile(a,[1,99.8]);return np.uint8(np.clip((a-lo)/max(hi-lo,1e-9)*255,0,255))
def extract(a):
 sift=cv2.SIFT_create(nfeatures=35000,contrastThreshold=.01,edgeThreshold=15)
 a=cv2.createCLAHE(clipLimit=2,tileGridSize=(16,16)).apply(a)
 return sift.detectAndCompute(a,None)
for sample in ['GSM9509172','GSM9509160','GSM9509170','GSM9509165']:
 t=time.time();p=next((OUT/'cache').glob(sample+'*ome.tif'));hp=next((OUT/'raw').glob('*hne_breast*'+('B4_B4' if sample=='GSM9509165' else 'B2_B2')+'*.ndpi'))
 with tifffile.TiffFile(p) as f:
  s=f.series[0];d=s.levels[3].asarray();ds=np.array(s.shape[::-1])/np.array(d.shape[::-1]);dx=ds*.2125
 d=cv2.resize(d,None,fx=dx[0]/2,fy=dx[1]/2,interpolation=cv2.INTER_AREA);d=norm(d);dk,dd=extract(d)
 he=openslide.OpenSlide(str(hp));mpp=float(he.properties['openslide.mpp-x']);lev=2;rgb=np.array(he.read_region((0,0),lev,he.level_dimensions[lev]).convert('RGB'));rgb=cv2.resize(rgb,None,fx=mpp*he.level_downsamples[lev]/2,fy=mpp*he.level_downsamples[lev]/2,interpolation=cv2.INTER_AREA)
 h=norm(rgb2hed(rgb)[:,:,0]);hk,hd=extract(h);pairs=cv2.BFMatcher().knnMatch(dd,hd,k=2);good=[m for m,n in pairs if m.distance<.78*n.distance]
 src=np.float32([dk[m.queryIdx].pt for m in good]);dst=np.float32([hk[m.trainIdx].pt for m in good]);matrix,inliers=cv2.estimateAffinePartial2D(src,dst,method=cv2.RANSAC,ransacReprojThreshold=8,maxIters=20000,confidence=.999,refineIters=20)
 if matrix is None:print(sample,'FAILED no matrix',flush=True);continue
 inliers=inliers.ravel().astype(bool);error=np.linalg.norm(src@matrix[:,:2].T+matrix[:,2]-dst,axis=1)*2
 print(sample,'matches',len(good),'inliers',inliers.sum(),'error_um',np.median(error[inliers]),'matrix',matrix.tolist(),'seconds',time.time()-t,flush=True)
 warp=cv2.warpAffine(d,matrix,(h.shape[1],h.shape[0]));mix=np.zeros((*h.shape,3),np.uint8);mix[:,:,0]=h;mix[:,:,1]=warp;mix[:,:,2]=h//2
 im=Image.fromarray(mix);im.thumbnail((2200,1600));im.save(DEST/(sample+'_overlay.png'));Image.fromarray(rgb).resize(im.size).save(DEST/(sample+'_HE.png'))
 # Conversion: DAPI micrometres -> H&E native pixels. Both registration arrays are 2 um/pixel.
 m=matrix.copy();m[:,:2]/=mpp;m[:,2]*=2/mpp
 rec={'sample':sample,'he_path':str(hp.relative_to(ROOT)),'dapi_path':str(p.relative_to(ROOT)),'he_mpp':mpp,'dapi_mpp':.2125,'registration_mpp':2,'matches':len(good),'inliers':int(inliers.sum()),'median_inlier_error_um':float(np.median(error[inliers])),'p95_inlier_error_um':float(np.quantile(error[inliers],.95)),'dapi_um_to_he_pixel':m.tolist(),'scale':float(np.linalg.norm(matrix[:,0])),'seconds':time.time()-t,'expression_accessed':False}
 (DEST/(sample+'.json')).write_text(json.dumps(rec,indent=2)+'\n');np.savez_compressed(DEST/(sample+'_landmarks.npz'),dapi_um=src*2,he_um=dst*2,inliers=inliers)
