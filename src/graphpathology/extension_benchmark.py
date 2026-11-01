"""Standardized CPU forward benchmark; run after concurrent experimental jobs finish."""
import json,time
import numpy as np
import pandas as pd
import torch
from .assets import sha
from .extension_evaluate import model_inventory
from .extension_models import ExtensionModel,graph_spec,model_inputs
from .extension_train import dataset,metadata
from .core import PeakRSS

def run(root):
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True);d=dataset(root,['TENX95'])['TENX95'];rows=[]
 for item in model_inventory(root):
  if item['origin']!='new':continue
  state=torch.load(root/item['checkpoint'],map_location='cpu',weights_only=True);model=ExtensionModel(item['variant']);model.load_state_dict(state['model']);model.eval();start=time.perf_counter();spec=graph_spec(d['xy_um'],'TENX95',item['seed'],item['variant']);graph_seconds=time.perf_counter()-start
  def forward():
   with torch.inference_mode():return model(torch.tensor(model_inputs(d['x'],spec,item['variant'])),spec)
  forward();forward();times=[]
  with PeakRSS() as rss:
   for _ in range(5):
    start=time.perf_counter();p=forward();times.append(time.perf_counter()-start)
  assert torch.isfinite(p).all()
  rows.append({'variant':item['variant'],'seed':item['seed'],'parameters':sum(p.numel() for p in model.parameters()),'n_locations':len(d['x']),'median_forward_seconds':float(np.median(times)),'min_forward_seconds':min(times),'max_forward_seconds':max(times),'graph_construction_seconds':graph_seconds,'peak_process_rss_mb':rss.peak/2**20,'threads':4,'sample':'TENX95','warmup':2,'repeats':5})
 out=root/'experiments/journal_revision/transformer_egnn/diagnostics';pd.DataFrame(rows).to_csv(out/'CPU_BENCHMARK.csv',index=False)
 report=metadata(root,status='COMPLETE',timing='Median5 full-section CPU forwards after2 warmups; includes model-specific feature aggregation/Tensor conversion, excludes graph construction and frozen CNN; process RSS includes interpreter/data; standalone after experimental workers finish',output_sha256=sha(out/'CPU_BENCHMARK.csv'))
 (out/'CPU_BENCHMARK.json').write_text(json.dumps(report,indent=2)+'\n');return report
