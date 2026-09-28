"""Additional graph diagnostics for the S2-280 principal geometry."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import sparse
from graph_controls import edge_array,lengths
from figures import LABELS,COLORS,save
fig,axes=plt.subplots(2,4,figsize=(17,7))
for j,sample in enumerate(IDS):
 xy=load(sample)['xy']
 for control in ['A3','A4','A5','A6']:
  matches=sorted((OUT/'graphs').glob(f'{sample}_knn6_cap250_{control}_*.json'));r=json.loads(matches[0].read_text());a=sparse.load_npz(ROOT/r['path']);le=np.sort(lengths(edge_array(a),xy));deg=np.diff(a.indptr)
  axes[0,j].hist(deg,bins=np.arange(0,17)-.5,density=True,histtype='step',label=LABELS[control],color=COLORS[control])
  axes[1,j].plot(le,np.arange(1,len(le)+1)/len(le),label=LABELS[control],color=COLORS[control])
 axes[0,j].set_title(sample);axes[0,j].set_xlabel('Node degree');axes[1,j].set_xlabel('Physical edge length (µm)');axes[1,j].set_xscale('log');axes[1,j].set_ylabel('Empirical CDF')
axes[0,0].legend(fontsize=8);fig.suptitle('S2-280 capped kNN6: first randomization; all four degree distributions overlap exactly');fig.tight_layout();save(fig,'graph_diagnostics_knn6.png')
