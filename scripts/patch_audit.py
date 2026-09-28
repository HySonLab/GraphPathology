"""Visual spot-check of cached H&E patches against coordinate-centered thumbnails."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'data/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import pandas as pd
from PIL import Image
from graphpathology.data import read_slide

def main():
    config=json.loads((ROOT/'configs/c0.json').read_text())
    meta=pd.read_csv(ROOT/'data/raw/HEST_v1_1_0.csv').set_index('id')
    ids=sum(config['splits'].values(),[])
    fig,axes=plt.subplots(4,6,figsize=(14,10))
    for row,sample in enumerate(ids):
        s=read_slide(ROOT/'data/raw',sample,float(meta.loc[sample,'pixel_size_um_estimated']))
        with h5py.File(s.patch_path) as f:
            width=float(f['img'].attrs['patch_size'])
            for j,fraction in enumerate([.2,.5,.8]):
                index=int(len(s.barcodes)*fraction)
                x,y=s.xy_px[index]*s.thumbnail_scale
                half=width*s.thumbnail_scale/2
                crop=Image.fromarray(s.thumbnail).crop((x-half,y-half,x+half,y+half))
                axes[row,2*j].imshow(crop)
                axes[row,2*j].set_title(f'{sample} {s.barcodes[index]}\nthumbnail crop',fontsize=9)
                axes[row,2*j+1].imshow(f['img'][index])
                axes[row,2*j+1].set_title('cached 224px patch',fontsize=9)
                axes[row,2*j].axis('off');axes[row,2*j+1].axis('off')
    fig.suptitle('Deterministic patch alignment spot-check: 20%, 50%, 80% of patch order\nThumbnail crops are low resolution; only gross tissue/color correspondence can be assessed.')
    fig.tight_layout(rect=(0,0,1,.94))
    fig.savefig(ROOT/'reports/C0/alignment/patch_spotcheck.png',dpi=150)
if __name__=='__main__': main()
