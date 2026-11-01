"""Pre-evaluation image-only affine intensity refinement diagnostic."""
exec(open('C7/refine_registration.py').read().split("for sample in ['GSM")[0])
for sample in ['GSM9509172','GSM9509160','GSM9509170','GSM9509165']:
 r=json.loads((dest/(sample+'_shape.json')).read_text());f=np.load(OUT/'cache'/(sample+'_registration_arrays.npz'));d=f['dapi'];h=f['he'];mat=np.array(r['matrix_2um']);w=cv2.warpAffine(d,mat,h.shape[::-1]);hm=cv2.resize(f['mask_he'],h.shape[::-1],interpolation=cv2.INTER_NEAREST);ys,xs=np.where(hm);x0=max(0,xs.min()-50);x1=min(h.shape[1],xs.max()+50);y0=max(0,ys.min()-50);y1=min(h.shape[0],ys.max()+50);a=h[y0:y1,x0:x1].astype('float32')/255;b=w[y0:y1,x0:x1].astype('float32')/255
 warp=np.eye(2,3,dtype='float32');history=[]
 for blur in [31,15,7,3]:
  try:cc,warp=cv2.findTransformECC(a,b,warp,cv2.MOTION_AFFINE,(cv2.TERM_CRITERIA_EPS|cv2.TERM_CRITERIA_COUNT,300,1e-6),None,blur);history.append({'blur':blur,'correlation':float(cc)})
  except cv2.error as e:history.append({'blur':blur,'error':str(e)[-150:]});break
 inv=np.linalg.inv(np.vstack([warp,[0,0,1]]));to=np.array([[1,0,x0],[0,1,y0],[0,0,1]]);final=(to@inv@np.linalg.inv(to)@np.vstack([mat,[0,0,1]]))[:2];w2=cv2.warpAffine(d,final,h.shape[::-1]);mix=np.stack([h,w2,h//2],-1);im=Image.fromarray(mix[y0:y1,x0:x1]);im.thumbnail((1600,1600));im.save(dest/(sample+'_ecc_overlay.png'));m=final.copy();m[:,:2]/=r['he_mpp'];m[:,2]*=2/r['he_mpp'];r.update(dapi_um_to_he_pixel=m.tolist(),matrix_2um=final.tolist(),ecc_history=history,status='DIAGNOSTIC_ONLY_NOT_QUALIFIED');(dest/(sample+'_ecc.json')).write_text(json.dumps(r,indent=2)+'\n');print(sample,history,flush=True)
