"""Qualify HER2ST using header/identifier metadata, supplied coordinates and H&E only."""
from pathlib import Path
import sys,json,gzip,io,hashlib
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C7';sys.path.insert(0,str(OUT/'deps'))
import pyzipper,numpy as np,pandas as pd
from PIL import Image
SAMPLES=[s+'1' for s in 'ABCDEFGH'];PWD=b'zNLXkYk3Q9znUseS';SPWD=b'yUx44SzG6NdB32gY'
PANELS=json.loads((ROOT/'C1/SELECTION_LOCK.json').read_text())['panels'];CACHE=OUT/'cache/her2';CACHE.mkdir(parents=True,exist_ok=True)
def count_metadata(z,s):
 name='count-matrices/'+s+'.tsv.gz'
 with gzip.GzipFile(fileobj=z.open(name,pwd=PWD)) as f:
  genes=f.readline().decode().rstrip('\r\n').split('\t')[1:]
  barcodes=[line.split(b'\t',1)[0].decode() for line in f if line.strip()]
 return genes,barcodes
if __name__=='__main__':
 assert not (OUT/'C7A/C7A_SPEC_LOCK.json').exists(),'C7A frozen'
 rows=[];genes_by={}
 with pyzipper.AESZipFile(OUT/'raw/count-matrices.zip') as cz,pyzipper.AESZipFile(OUT/'raw/spot-selections.zip') as sz:
  for s in SAMPLES:
   genes,bc=count_metadata(cz,s);genes_by[s]=genes
   df=pd.read_csv(gzip.GzipFile(fileobj=sz.open(s+'_selection.tsv.gz',pwd=SPWD)),sep='\t');df['barcode']=df.x.astype(int).astype(str)+'x'+df.y.astype(int).astype(str);df=df[df.selected==1].set_index('barcode');available=[b for b in bc if b in df.index];df=df.loc[available];assert len(available)==len(bc) and len(set(bc))==len(bc)
   a=np.column_stack([df.new_x,df.new_y,np.ones(len(df))]);coef=np.linalg.lstsq(a,df[['pixel_x','pixel_y']],rcond=None)[0];singular=np.linalg.svd(coef[:2],compute_uv=False);mpp=float(200/np.mean(singular));res=float(np.max(np.abs(a@coef-df[['pixel_x','pixel_y']].values)));assert res<2 and singular.max()/singular.min()<1.05
   np.savez_compressed(CACHE/(s+'_metadata.npz'),barcodes=np.array(bc),xy_px=df[['pixel_x','pixel_y']].values,xy_um=df[['pixel_x','pixel_y']].values*mpp,mpp=mpp)
   rows.append({'sample_id':s,'patient_id':'HER2ST_'+s[0],'num_locations':len(bc),'num_assay_genes':len(genes),'mpp':mpp,'coordinate_fit_max_error_px':res,'coordinate_anisotropy':float(singular.max()/singular.min()),'overlaps':{p:len(set(v)&set(genes)) for p,v in PANELS.items() if p!='S1'}})
 common=set.intersection(*(set(g) for g in genes_by.values()));retained={p:[g for g in PANELS[p] if g in common] for p in ['S2-50','S2-280']};assert len(retained['S2-50'])>=40
 record={'samples':rows,'panel_intersections':retained,'gene_headers':genes_by,'minimum_overlap_S2_50':.8,'expression_values_loaded':False,'source_filters':'Published processed gene headers and selected spot identities; no external expression values parsed for qualification.'};(OUT/'C7A/logs/HER2ST_METADATA_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(rows,indent=2))
 if (OUT/'raw/images.zip').exists():
  with pyzipper.AESZipFile(OUT/'raw/images.zip') as z:
   print('image paths',z.namelist()[:5])
   for s in SAMPLES:
    n=next(n for n in z.namelist() if n.endswith('/HE/'+s+'.jpg') or n=='HE/'+s+'.jpg');p=CACHE/(s+'.jpg')
    if not p.exists():p.write_bytes(z.read(n,pwd=PWD))
    print(s,Image.open(p).size)
