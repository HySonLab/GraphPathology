"""Metadata-only image registration contact sheet before expression evaluation."""
from contract10 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
fig,axs=plt.subplots(4,6,figsize=(18,12))
for ax,s in zip(axs.flat,SAMPLES):
 im=Image.open(ROOT/cfg()['inputs'][s]['image'])
 with np.load(CACHE/(s+'_metadata.npz')) as f:xy=f['xy_px']*float(f['hires_scale'])
 ax.imshow(im);ax.scatter(xy[:,0],xy[:,1],s=.3,c='cyan',alpha=.35);ax.set_title(s,fontsize=10);ax.axis('off')
for ax in list(axs.flat)[len(SAMPLES):]:ax.axis('off')
fig.suptitle('C10 metadata-only coordinate overlays: first section of each patient')
fig.tight_layout();fig.savefig(A/'logs/COORDINATE_QA.png',dpi=130);plt.close(fig)
