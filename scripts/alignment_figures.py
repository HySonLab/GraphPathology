"""Show barcode-joined ST positions against embedded H&E thumbnails."""
import json
from pathlib import Path
import sys
import os
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[1]/'data/matplotlib'))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from graphpathology.data import read_slide
ROOT = Path(__file__).resolve().parents[1]

def main():
    config = json.loads((ROOT/'configs/c0.json').read_text())
    meta = pd.read_csv(ROOT/'data/raw/HEST_v1_1_0.csv').set_index('id')
    samples = sum(config['splits'].values(), [])
    fig, axes = plt.subplots(2,4,figsize=(16,10))
    for col, sample in enumerate(samples):
        s = read_slide(ROOT/'data/raw',sample,float(meta.loc[sample,'pixel_size_um_estimated']))
        for row in range(2):
            ax = axes[row,col]
            ax.imshow(s.thumbnail)
            if row:
                ax.scatter(*(s.xy_px*s.thumbnail_scale).T,s=.3,c='#00ffff',alpha=.35)
            ax.set_title(f'{sample}: '+ ('joined spots' if row else 'H&E'))
            ax.axis('off')
    fig.suptitle('C0 alignment audit: HEST IDC\nCoordinates joined by barcode; patch row/column converted to x/y. Thumbnail view cannot certify cell-level registration.')
    fig.tight_layout()
    out=ROOT/'reports/C0/alignment'; out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/'cohort_overlay.png',dpi=160)
    plt.close(fig)
if __name__=='__main__': main()
