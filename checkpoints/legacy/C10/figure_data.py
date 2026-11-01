"""Export frozen spatial maps and label-free schematic source coordinates."""
from synthesis import *
def main():
 d=read('C5/C5_PANEL_METRICS.csv');d=d[(d.phase=='primary')&(d.panel=='S2-50')&(d.nominal_coverage==.2)&(d.replicate==0)&(d.model=='B3')&(d.encoder==ENC)&d.control.isna()];rows=[]
 with np.load(ROOT/'data/cache/TENX99_targets.npz') as f:y=f['normalized_targets'];bc=f['barcodes'];xy=f['xy_um'];gg=f['genes']
 for _,r in d.iterrows():
  with np.load(ROOT/r.prediction_path) as f:
   ix=pd.Index(bc).get_indexer(f['hidden_barcodes']);g=str(f['genes'][0]);j=list(gg).index(g);err=abs(y[ix,j]-f['prediction'][:,0]);ob=pd.Index(bc).get_indexer(f['observed_barcodes'])
  for k,e in zip(ix,err):rows.append({'regime':r.regime,'x_um':xy[k,0],'y_um':xy[k,1],'absolute_error':e,'observed':False,'gene':g,'source_file':r.prediction_path,'source_sha256':r.prediction_sha256,'source_barcode':bc[k],'source_gene_index':0})
  for k in ob:rows.append({'regime':r.regime,'x_um':xy[k,0],'y_um':xy[k,1],'absolute_error':np.nan,'observed':True,'gene':g,'source_file':r.prediction_path,'source_sha256':r.prediction_sha256,'source_barcode':bc[k],'source_gene_index':0})
 save(pd.DataFrame(rows),'Fig3_map_source')
 xy=np.array([(x,y) for x in range(12) for y in range(10)],float);out=[]
 for regime in ['V1','V2']:
  forbidden=(np.linalg.norm(xy-[6,5],axis=1)<3) if regime=='V2' else np.zeros(len(xy),bool);ids=np.random.default_rng(17).permutation(np.flatnonzero(~forbidden))[:24]
  for i,(x,y) in enumerate(xy):out.append({'regime':regime,'x':x,'y':y,'observed':i in ids,'forbidden':forbidden[i],'source':'Synthetic illustrative coordinates; not experimental results; seed17'})
 save(pd.DataFrame(out),'Fig1_source')
if __name__=='__main__':main()
