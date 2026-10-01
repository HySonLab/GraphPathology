"""Visual QA of supplied HER2ST spot-to-H&E coordinates, before outcomes."""
from support import *
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
Image.MAX_IMAGE_PIXELS=120_000_000
fig,axes=plt.subplots(2,4,figsize=(16,9))
for ax,s in zip(axes.ravel(),SAMPLES):
 im=Image.open(CACHE/(s+'.jpg'));size=im.size;im.thumbnail((1100,1100));ax.imshow(im,extent=[0,size[0],size[1],0])
 with np.load(CACHE/(s+'_metadata.npz')) as f:xy=f['xy_px'];mpp=float(f['mpp'])
 ax.scatter(xy[:,0],xy[:,1],s=3,facecolors='none',edgecolors='lime',linewidth=.4);ax.set_title(f'{s}: {len(xy)} spots; {mpp:.3f} µm/px');ax.set_axis_off()
fig.suptitle('External HER2ST: supplied H&E coordinates (all selected spots)\nFirst section of each independent patient; no expression values used');fig.tight_layout();fig.savefig(A/'logs/HER2ST_ALIGNMENT_QA.png',dpi=150);plt.close(fig)
