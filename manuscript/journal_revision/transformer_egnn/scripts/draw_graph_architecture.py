#!/usr/bin/env python3
"""Rebuild Figure 2 using synthetic spatial coordinates; no research data or training used."""
from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay, cKDTree

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, default=root/'figures',
                    help='Directory for Figure2.pdf and Figure2.svg.')
parser.add_argument('--data-dir', type=Path, default=root/'source_data',
                    help='Directory for synthetic coordinates and edge CSVs.')
args = parser.parse_args()
figdir = args.output_dir
figdir.mkdir(parents=True, exist_ok=True)
args.data_dir.mkdir(parents=True, exist_ok=True)
rng=np.random.default_rng(23)
xx,yy=np.meshgrid(np.linspace(0,800,6),np.linspace(0,800,6))
p=np.stack([xx.flatten(),yy.flatten()],axis=-1)+rng.uniform(-45,45,(36,2))
n=len(p); dm=np.linalg.norm(p[:,None,:]-p[None,:,:],axis=-1)

def make_graph(kind):
    e=set()
    if kind=='kNN (k=6, capped at 250 µm)':
        ind=np.argsort(dm,axis=1)[:,1:7]
        for i in range(n):
            for j in ind[i]:
                if dm[i,j] <= 250:
                    e.add(tuple(sorted((i,int(j)))))
    elif kind.startswith('Radius'):
        cut=150 if '150' in kind else 250
        e={(i,j) for i in range(n) for j in range(i+1,n) if dm[i,j]<=cut}
    else:
        tris=Delaunay(p).simplices
        pairs={tuple(sorted((int(tri[u]),int(tri[(u+1)%3])))) for tri in tris for u in range(3)}
        e={edge for edge in pairs if dm[edge[0],edge[1]]<=250}
    return sorted(e)

labels=['kNN (k=6, capped at 250 µm)','Radius (150 µm)','Radius (250 µm)','Delaunay (capped at 250 µm)']
fig=plt.figure(figsize=(11.8,10.2)); gs=fig.add_gridspec(3,2,height_ratios=[1,1,0.59],hspace=0.40,wspace=0.20)
for i,label in enumerate(labels):
    ax=fig.add_subplot(gs[i//2,i%2]); edges=make_graph(label)
    for a,b in edges:
        ax.plot([p[a,0],p[b,0]],[p[a,1],p[b,1]],color='#86a8b6',alpha=.62,lw=.95,zorder=1)
    ax.scatter(p[:,0],p[:,1],s=36,color='#176d83',edgecolors='white',linewidths=.5,zorder=3)
    ax.set_aspect('equal');ax.set_xlim(-100,900);ax.set_ylim(-100,900)
    ax.set_title(f'{"abcd"[i]}. {label}  |  {len(edges)} edges',size=11,loc='left',pad=9)
    ax.set_xlabel('Physical x (µm)',size=9);ax.set_ylabel('Physical y (µm)',size=9)
    ax.grid(False);ax.tick_params(labelsize=8)
ax=fig.add_subplot(gs[2,:]);ax.axis('off')
steps=[(0.005,0.17,'Tissue locations\n+ physical coordinates'),(.265,.17,'Coordinate-defined\nadjacency / fixed operator'),(.525,.17,'Image embedding\n+ spatial aggregation'),(.785,.17,'Trained regression head\n→ 50-gene predictions')]
for x,y,s in steps:
    ax.add_patch(plt.Rectangle((x,y),.205,.56,transform=ax.transAxes,facecolor='#e8f3f4',edgecolor='#35788a',linewidth=1.2))
    ax.text(x+.1025,y+.28,s,transform=ax.transAxes,ha='center',va='center',fontsize=10)
for i in range(3):
    ax.annotate('',xy=(steps[i+1][0]-.008,.45),xytext=(steps[i][0]+.213,.45),xycoords='axes fraction',textcoords='axes fraction',arrowprops={'arrowstyle':'-|>','lw':1.7,'color':'#406b72'})
ax.text(.01,.91,'e. Geometry-aware contextual representation pathway (development models)',transform=ax.transAxes,fontsize=12,weight='bold')
ax.text(.005,.04,'Illustrative coordinates only. The evaluation-stage molecular-anchor graph is constructed separately from observed sites.',transform=ax.transAxes,fontsize=9)
fig.savefig(figdir/'Figure2.pdf',bbox_inches='tight');fig.savefig(figdir/'Figure2.svg',bbox_inches='tight');
with open(args.data_dir/'Figure2_synthetic_locations.csv','w') as f:
    f.write('location_id,x_um,y_um\n')
    for i,(x,y) in enumerate(p): f.write(f'{i},{x:.6f},{y:.6f}\n')
with open(args.data_dir/'Figure2_synthetic_edges.csv','w') as f:
    f.write('graph,source_id,target_id,distance_um\n')
    for lab in labels:
        for a,b in make_graph(lab): f.write(f'"{lab}",{a},{b},{dm[a,b]:.6f}\n')
plt.close(fig)
print('Wrote Figure2.pdf/SVG and synthetic source CSVs')
