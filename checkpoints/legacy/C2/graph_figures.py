"""Coordinate-only visual audit of the four frozen graph candidates."""
import os,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C2'
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'data/matplotlib'));os.environ.setdefault('XDG_CACHE_HOME',str(ROOT/'data/.cache'))
import numpy as np
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
with np.load(ROOT/'data/cache/NCBI783_targets.npz') as f:xy=f['xy_um']
center=np.median(xy,axis=0);keep=np.max(abs(xy-center),axis=1)<=650
fig,axes=plt.subplots(1,4,figsize=(14,4))
config=json.loads((OUT/'configs/protocol.json').read_text())
for ax,g in zip(axes,config['graphs']):
 a=sparse.load_npz(OUT/'graphs'/f'NCBI783_{g["name"]}.npz');i,j=a.nonzero();valid=keep[i]&keep[j]&(i<j)
 segments=np.stack([xy[i[valid]],xy[j[valid]]],axis=1)/1000
 ax.add_collection(LineCollection(segments,colors='#94a3b8',linewidths=.6));ax.scatter(*xy[keep].T/1000,s=9,c='#0369a1')
 ax.set_xlim((center[0]-650)/1000,(center[0]+650)/1000);ax.set_ylim((center[1]+650)/1000,(center[1]-650)/1000);ax.set_aspect('equal');ax.set(title=g['name'],xlabel='x (mm)',ylabel='y (mm)')
fig.suptitle('Graph audit: fixed coordinate-centered training-slide window (NCBI783)\nEdges shown only when both endpoints lie in the display window; full graphs are slide-local.')
fig.tight_layout(rect=(0,0,1,.85));fig.savefig(OUT/'figures/graph_construction_audit.png',dpi=160,bbox_inches='tight')
