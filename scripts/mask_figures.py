"""Coordinate-only diagnostic figures; no expression labels are read for masks."""
import json
import os
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'data/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from graphpathology.data import read_slide
from graphpathology.masks import make_mask,stable_seed,REGIMES

def main():
    config=json.loads((ROOT/'configs/c0.json').read_text())
    meta=pd.read_csv(ROOT/'data/raw/HEST_v1_1_0.csv').set_index('id')
    ids=sum(config['splits'].values(),[])
    fig,axes=plt.subplots(len(ids),3,figsize=(12,14))
    for row,sample in enumerate(ids):
        slide=read_slide(ROOT/'data/raw',sample,float(meta.loc[sample,'pixel_size_um_estimated']))
        for col,regime in enumerate(REGIMES):
            mask,info=make_mask(slide.xy_um,regime,.2,stable_seed(sample,regime,.2,0,config['seed']))
            xy=slide.xy_um/1000
            ax=axes[row,col]
            ax.scatter(*xy[~mask].T,s=1,c='#d1d5db',label='Hidden')
            ax.scatter(*xy[mask].T,s=1,c='#075985',label='Observed')
            ax.set_title(f'{sample}: {regime}\nobserved {mask.mean():.1%}')
            ax.set_aspect('equal'); ax.invert_yaxis(); ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    fig.suptitle('Frozen S3 mask examples at nominal 20% spot coverage\nSquare block is hidden; circular capture window is observed. Tissue-area coverage is not estimated.')
    fig.tight_layout(rect=(0,0,1,.96))
    fig.savefig(ROOT/'reports/C0/alignment/mask_geometries.png',dpi=130)
    plt.close(fig)
if __name__=='__main__': main()
