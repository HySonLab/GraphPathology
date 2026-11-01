import json
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from graphpathology.backbone_data import record,write
from graphpathology.assets import sha
root=Path(__file__).resolve().parents[1];base=root/'experiments/journal_revision/backbone_graph_comparison';rows=json.loads((base/'ALL_CANDIDATE_FIT_MANIFEST.json').read_text());selected=json.loads((base/'SELECTED_FIT_MANIFEST.json').read_text());out=base/'histories';out.mkdir(exist_ok=True)
report=[]
for r in rows:
    h=pd.DataFrame(r['history']) if 'history' in r else pd.read_csv(root/Path(r['source_fit_record']).with_name(Path(r['source_fit_record']).stem+'_history.csv'))
    assert len(h)==r['epochs_ran'] and np.isfinite(h[['training_mse','validation_mse']].to_numpy()).all()
    assert abs(float(h.iloc[r['selected_epoch']-1].validation_mse)-r['validation_mse'])<1e-12
    p=out/(r['run_id']+'_history.csv');h.to_csv(p,index=False)
    report.append({'run_id':r['run_id'],'history':str(p.relative_to(root)),'history_sha256':sha(p),'epochs_ran':len(h),'selected_epoch':r['selected_epoch'],'at_budget_ceiling':len(h)==80,'all_finite':True})
fig,axes=plt.subplots(2,6,figsize=(18,6),squeeze=False)
for i,backbone in enumerate(['resnet18','phikon']):
 for j,(model,control) in enumerate([('local','local'),('MEAN','spatial'),('GCN','spatial'),('SAGE','spatial'),('GCN','shuffled'),('SAGE','shuffled')]):
  ax=axes[i,j]
  for r in [r for r in selected if r['backbone']==backbone and r['model']==model and r['control']==control]:
   h=pd.read_csv(out/(r['run_id']+'_history.csv'));line,=ax.plot(h.epoch,h.validation_mse,label=f"seed {r['seed']}");ax.plot(h.epoch,h.training_mse,':',color=line.get_color(),alpha=.5);ax.plot(r['selected_epoch'],r['validation_mse'],'o',color=line.get_color(),ms=3)
  ax.set_title(f'{backbone}: {model} {control}',fontsize=9);ax.set_xlabel('Epoch');ax.set_ylabel('MSE');ax.legend(fontsize=6);ax.grid(alpha=.2)
fig.suptitle('Validation-selected heads: solid validation, dotted source fitting; circles retained epochs',fontsize=12);fig.tight_layout();(base/'diagnostics').mkdir(exist_ok=True);fig.savefig(base/'diagnostics/SelectedLearningCurves.png',dpi=130);plt.close(fig)
write(base/'SELECTED_LEARNING_CURVE_CHECKS.json',record(root,status='COMPLETE',candidate_histories_checked=72,selected_heads=36,finite_losses=True,retained_validation_loss_matches=True,all_histories_preserved=True,records=report))
