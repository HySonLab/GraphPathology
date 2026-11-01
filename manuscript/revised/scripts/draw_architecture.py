#!/usr/bin/env python3
"""Draw the verified neural architecture and study schematic.

Architecture and training settings were checked against the public analysis
revision 03506defbe205cabd13865780d6a1224b348a781. The source map is in
METHODS_SOURCE_AUDIT.txt. Spatial examples reuse the supplied synthetic data.
No fitted models or scientific result arrays are modified by this script.
"""
from pathlib import Path
import csv
import io
import os
import tempfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9,
                     'pdf.fonttype':42, 'svg.fonttype':'none'})
fig = plt.figure(figsize=(7.2, 7.6), facecolor='white')
ax = fig.add_axes([0,0,1,1]); ax.set(xlim=(0,7.2), ylim=(0,7.6)); ax.axis('off')
ink='#1d2935'; teal='#087e8b'; rust='#a54b2a'; violet='#8050a0'

def text(x,y,s,size=9,ha='left',va='center',weight='normal',color=ink,**kw):
    return ax.text(x,y,s,fontsize=size,ha=ha,va=va,fontweight=weight,color=color,**kw)

def box(x,y,w,h,s,fc='#f1f5f7',ec='#8d9ba6',size=9,lw=.9):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.025,rounding_size=0.06',
                               facecolor=fc,edgecolor=ec,linewidth=lw))
    text(x+w/2,y+h/2,s,size,ha='center',linespacing=1.35)

def arrow(x1,y1,x2,y2,color=ink):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',
                 mutation_scale=9,linewidth=1,color=color,shrinkA=1,shrinkB=1))

def heading(y,letter,title):
    text(.1,y,letter,12,weight='bold')
    text(.35,y,title,10.1,weight='bold')

heading(7.40,'a','Cohort roles and observation geometry')
rows=list(csv.DictReader((ROOT/'source_data/Fig1_source.csv').open()))
for left,reg,title in [(.16,'V1','Distributed anchors (V1)'),(2.05,'V2','Contiguous region (V2)')]:
    pa=fig.add_axes([left/7.2,5.70/7.6,1.68/7.2,1.35/7.6])
    selected=[r for r in rows if r['regime']==reg]
    colors=[teal if r['observed'].lower()=='true' else
            violet if r['forbidden'].lower()=='true' else '#cbd2d9' for r in selected]
    pa.scatter([float(r['x']) for r in selected],[float(r['y']) for r in selected],
               s=13,c=colors,linewidths=0)
    pa.set_aspect('equal');pa.set_xticks([]);pa.set_yticks([])
    for spine in pa.spines.values():spine.set_visible(False)
    pa.set_title(title,fontsize=8.8,pad=4,color=ink)

text(4.04,6.94,'Development: HEST, 4 patients',9.1,weight='bold')
text(4.04,6.66,'2 train / 1 validation / 1 evaluation',8.4)
text(4.04,6.32,'Transportability: HER2ST, 8 patients',9.1,weight='bold')
text(4.04,6.04,'Reused for calibration development',8.4)
text(4.04,5.70,'Fresh evaluation: Wu 3 + Bassiouni 22',9.1,weight='bold')
text(4.04,5.43,'Three separate Wu patients supply calibration',8.3)
for x,c,label in [(.25,teal,'Observed'),(1.43,'#cbd2d9','Hidden'),(2.47,violet,'Scored region')]:
    ax.scatter([x],[5.43],s=24,c=c,linewidths=0)
    text(x+.13,5.43,label,8.5)
ax.plot([.1,7.1],[5.18,5.18],color='#cad3da',lw=.8)

heading(4.97,'b','Frozen encoder and three trained regression heads')
box(.12,3.79,1.24,.90,'Nominal 112 μm\nH&E patch\n224 × 224 RGB',size=8.6)
box(1.57,3.69,1.95,1.10,'Frozen ResNet18\n4 stages: 64/128/256/512\n2 basic blocks per stage\nGlobal average pooling\n512-dimensional vector',
    fc='#eef5f8',ec='#688ba5',size=8.2)
box(3.73,3.79,1.18,.90,'Feature\nstandardization\n(fit on training)',fc='#eef5f8',ec='#688ba5',size=8.2)
box(5.12,3.69,1.96,1.10,'Three MLP heads\nLinear: 512 → 256\nGELU → dropout 0.1\nLinear: 256 → 50\nSeeds: 0, 1, 2',
    fc='#f4f0f8',ec='#9b83ad',size=8.6)
arrow(1.38,4.24,1.55,4.24);arrow(3.54,4.24,3.71,4.24);arrow(4.93,4.24,5.10,4.24)
text(.12,3.40,'Patch processing: LANCZOS → 224; official resize → 256;\ncenter crop → 224; ImageNet channel normalization',8.0,va='center')
text(6.10,3.36,'Invert target scaling for each head;\nmean of three 50-gene predictions',8.0,ha='center')
arrow(6.10,3.68,6.10,3.57,color=violet)
ax.plot([.1,7.1],[3.05,3.05],color='#cad3da',lw=.8)

heading(2.83,'c','Molecular interpolation and fixed fusion')
box(.12,1.82,1.40,.69,'Coordinates +\nobserved expression',fc='#eaf4f3',ec=teal,size=8.6)
box(1.79,1.68,1.88,.97,'Up to 16 anchors\nV1: inverse distance²\nV2: equal weights',fc='#eaf4f3',ec=teal,size=8.8)
box(3.94,1.82,1.28,.69,'Molecular\nprediction',fc='#eaf4f3',ec=teal,size=9)
box(5.58,1.67,1.50,.99,'Fixed fusion\n0.8 molecular\n+ 0.2 image',fc='#faf0eb',ec=rust,size=9)
arrow(1.54,2.17,1.77,2.17,teal);arrow(3.69,2.17,3.92,2.17,teal);arrow(5.24,2.17,5.56,2.17,teal)
text(4.58,1.56,r'$\widehat{\mathbf{y}}^{\mathrm{ST}}_i$',10,ha='center',color=teal)
box(.12,.65,1.40,.66,'Image-only\nprediction',fc='#faf0eb',ec=rust,size=9)
box(1.79,.65,3.43,.66,'Frozen source-to-Wu image calibration\nTransported unchanged to Bassiouni',fc='#faf0eb',ec=rust,size=8.8)
arrow(1.54,.98,1.77,.98,rust)
ax.plot([5.24,6.33],[.98,.98],color=rust,lw=1)
arrow(6.33,.98,6.33,1.65,rust)
text(3.58,.25,'Score hidden locations only: image, molecular and fusion predictions',9,ha='center',weight='bold')
text(3.58,.03,'V1 scores all hidden locations; V2 scores the designated contiguous region.',8.3,ha='center')

out=ROOT/'figures';out.mkdir(exist_ok=True)
def atomic_save(path, format, **kwargs):
    buffer = io.BytesIO()
    fig.savefig(buffer, format=format, **kwargs)
    payload = buffer.getvalue()
    if format == 'pdf':
        assert payload.rstrip().endswith(b'%%EOF'), 'Incomplete PDF export'
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)

for suffix in ['pdf','svg']:
    atomic_save(out/f'Figure1.{suffix}', suffix, bbox_inches='tight', pad_inches=.02)
atomic_save(out/'Figure1_preview.png', 'png', dpi=180, bbox_inches='tight', pad_inches=.02)
plt.close(fig)
print('Saved Figure1.pdf, Figure1.svg and preview.')
