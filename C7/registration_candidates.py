"""Explore label-free contrast and orientation choices before qualification lock."""
exec(open('C7/register.py').read().split("for sample in ['GSM")[0])
for sample in ['GSM9509172','GSM9509160','GSM9509170','GSM9509165']:
 p=next((OUT/'cache').glob(sample+'*ome.tif'));hp=next((OUT/'raw').glob('*hne_breast*'+('B4_B4' if sample=='GSM9509165' else 'B2_B2')+'*.ndpi'))
 with tifffile.TiffFile(p) as f:
  s=f.series[0];d=s.levels[3].asarray();ds=np.array(s.shape[::-1])/np.array(d.shape[::-1]);dx=ds*.2125
 d=cv2.resize(d,None,fx=dx[0]/2,fy=dx[1]/2,interpolation=cv2.INTER_AREA);d=norm(d)
 he=openslide.OpenSlide(str(hp));mpp=float(he.properties['openslide.mpp-x']);rgb=np.array(he.read_region((0,0),2,he.level_dimensions[2]).convert('RGB'));rgb=cv2.resize(rgb,None,fx=mpp*4/2,fy=mpp*4/2,interpolation=cv2.INTER_AREA)
 hs=[norm(rgb2hed(rgb)[:,:,0]),norm(255-rgb[:,:,0]),norm(255-cv2.cvtColor(rgb,cv2.COLOR_RGB2GRAY))];best=None
 for contrast,h in enumerate(hs):
  hk,hd=extract(h)
  for flip in [False,True]:
   d0=np.fliplr(d).copy() if flip else d;dk,dd=extract(d0);pairs=cv2.BFMatcher().knnMatch(dd,hd,k=2);good=[m for m,n in pairs if m.distance<.88*n.distance]
   src=np.float32([dk[m.queryIdx].pt for m in good]);dst=np.float32([hk[m.trainIdx].pt for m in good]);mat,ins=cv2.estimateAffinePartial2D(src,dst,method=cv2.RANSAC,ransacReprojThreshold=6,maxIters=40000,confidence=.999,refineIters=20)
   if mat is None:continue
   ins=ins.ravel().astype(bool);scale=np.linalg.norm(mat[:,0]);score=int(ins.sum()) if .8<scale<1.2 else 0
   print(sample,contrast,flip,'inliers',ins.sum(),'scale',scale,flush=True)
   if best is None or score>best[0]:best=(score,contrast,flip,mat,src,dst,ins,h)
 score,contrast,flip,mat,src,dst,ins,h=best
 if flip:
  src[:,0]=d.shape[1]-1-src[:,0];mat=np.vstack([mat,[0,0,1]])@np.array([[-1,0,d.shape[1]-1],[0,1,0],[0,0,1]]);mat=mat[:2]
 error=np.linalg.norm(src@mat[:,:2].T+mat[:,2]-dst,axis=1)*2;warp=cv2.warpAffine(d,mat,(h.shape[1],h.shape[0]));mix=np.stack([h,warp,h//2],axis=-1);im=Image.fromarray(mix);im.thumbnail((2200,1600));im.save(DEST/(sample+'_candidate_overlay.png'))
 m=mat.copy();m[:,:2]/=mpp;m[:,2]*=2/mpp
 r={'sample':sample,'he_path':str(hp.relative_to(ROOT)),'he_mpp':mpp,'inliers':score,'median_inlier_error_um':float(np.median(error[ins])),'dapi_um_to_he_pixel':m.tolist(),'contrast':contrast,'flip':flip,'scale':float(np.linalg.norm(mat[:,0])),'expression_accessed':False}
 (DEST/(sample+'_candidate.json')).write_text(json.dumps(r,indent=2)+'\n');np.savez_compressed(DEST/(sample+'_candidate_landmarks.npz'),dapi_um=src*2,he_um=dst*2,inliers=ins)
