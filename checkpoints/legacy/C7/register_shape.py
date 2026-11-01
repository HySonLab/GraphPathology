"""Image-only tissue contour registration diagnostic, not an accepted alignment."""
exec(open('C7/register.py').read().split("for sample in ['GSM")[0])
from scipy.ndimage import binary_fill_holes,distance_transform_edt
from scipy.optimize import minimize

def mask(a):
 z=(a>25).astype('uint8');z=cv2.morphologyEx(z,cv2.MORPH_CLOSE,np.ones((15,15),np.uint8));n,l,st,c=cv2.connectedComponentsWithStats(z)
 return z
for sample in ['GSM9509172','GSM9509160','GSM9509170','GSM9509165']:
 p=next((OUT/'cache').glob(sample+'*ome.tif'));hp=next((OUT/'raw').glob('*hne_breast*'+('B4_B4' if sample=='GSM9509165' else 'B2_B2')+'*.ndpi'))
 with tifffile.TiffFile(p) as f:
  s=f.series[0];d=s.levels[3].asarray();ds=np.array(s.shape[::-1])/np.array(d.shape[::-1]);dx=ds*.2125
 d=cv2.resize(d,None,fx=dx[0]/2,fy=dx[1]/2,interpolation=cv2.INTER_AREA);d=norm(d)
 he=openslide.OpenSlide(str(hp));mpp=float(he.properties['openslide.mpp-x']);rgb=np.array(he.read_region((0,0),2,he.level_dimensions[2]).convert('RGB'));rgb=cv2.resize(rgb,None,fx=mpp*2,fy=mpp*2,interpolation=cv2.INTER_AREA);h=norm(rgb2hed(rgb)[:,:,0])
 # Target image region determined only by connected tissue components and gross morphology.
 hs=cv2.resize(rgb,None,fx=.2,fy=.2);sat=cv2.cvtColor(hs,cv2.COLOR_RGB2HSV)[:,:,1];hm=mask(sat);hm=binary_fill_holes(hm).astype('uint8');n,l,st,c=cv2.connectedComponentsWithStats(hm)
 ids=[i for i in range(1,n) if st[i,4]>3000];ids.sort(key=lambda i:c[i,0])
 if sample=='GSM9509165':ids=[i for i in ids if c[i,0]<hm.shape[1]*.55]
 else:ids=[ids[['GSM9509172','GSM9509160','GSM9509170'].index(sample)]]
 hm=np.uint8(np.isin(l,ids));dm=mask(cv2.resize(d,None,fx=.2,fy=.2));n,l,st,c=cv2.connectedComponentsWithStats(dm);ids=[i for i in range(1,n) if st[i,4]>max(st[1:,4])*.05];dm=binary_fill_holes(np.isin(l,ids)).astype('uint8')
 dp=np.argwhere(dm)[:,::-1];hp0=np.argwhere(hm)[:,::-1];dc=dp.mean(0);hc=hp0.mean(0);_,dv=np.linalg.eigh(np.cov(dp.T));_,hv=np.linalg.eigh(np.cov(hp0.T));da=np.arctan2(dv[1,-1],dv[0,-1]);ha=np.arctan2(hv[1,-1],hv[0,-1]);scale=np.sqrt(hm.sum()/dm.sum());best=None
 # Optimizing mask overlap is a label-free preliminary image registration only.
 for flip in [1,-1]:
  for k in [0,1]:
   theta=ha-(da if flip==1 else np.pi-da)+k*np.pi
   def matrix(v):
    a=v[0];s=np.exp(v[1]);r=s*np.array([[np.cos(a)*flip,-np.sin(a)],[np.sin(a)*flip,np.cos(a)]]);return np.column_stack([r,hc+v[2:]-r@dc])
   def fun(v):
    w=cv2.warpAffine(dm,matrix(v),hm.shape[::-1],flags=cv2.INTER_NEAREST);return -2.0*float(np.sum(w*hm))/(w.sum()+hm.sum()+1)
   r=minimize(fun,[theta,np.log(scale),0,0],method='Powell',options={'maxiter':120,'xtol':1e-4,'ftol':1e-5});score=-r.fun
   if best is None or score>best[0]:best=(score,matrix(r.x),flip)
 score,mat,flip=best;mat[:,2]*=5;w=cv2.warpAffine(d,mat,h.shape[::-1]);mix=np.stack([h,w,h//2],axis=-1);im=Image.fromarray(mix);im.thumbnail((2200,1600));im.save(DEST/(sample+'_shape_overlay.png'));m=mat.copy();m[:,:2]/=mpp;m[:,2]*=2/mpp
 r={'sample':sample,'he_path':str(hp.relative_to(ROOT)),'he_mpp':mpp,'mask_dice':score,'dapi_um_to_he_pixel':m.tolist(),'matrix_2um':mat.tolist(),'flip':flip,'scale':float(np.linalg.norm(mat[:,0])),'expression_accessed':False,'status':'DIAGNOSTIC_ONLY_NOT_QUALIFIED'};(DEST/(sample+'_shape.json')).write_text(json.dumps(r,indent=2)+'\n');np.savez_compressed(OUT/'cache'/(sample+'_registration_arrays.npz'),dapi=d,he=h,rgb=rgb,mask_he=hm,mask_dapi=dm);print(sample,r,flush=True)
