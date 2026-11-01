"""Read-only pre-test replay of completed validation-only model selection."""
from common import *
import torch
from threadpoolctl import threadpool_limits
spec=importlib.util.spec_from_file_location('c4_runtime',OUT/'run.py');run=importlib.util.module_from_spec(spec);spec.loader.exec_module(run)
def main():
 lock=json.loads((A/'SELECTION_LOCK.json').read_text());assert not (A/'logs/TEST_ACCESS.json').exists()
 for p,h in lock['source_sha256'].items():assert sha(ROOT/p)==h
 records=[json.loads(p.read_text()) for p in (A/'logs/runs').glob('*.json')];assert len(records)==116
 count=0
 for encoder in run.CFG['encoders']:
  data={s:run.load(encoder,s) for s in SPLITS['splits']['train']+SPLITS['splits']['validation']};sample=SPLITS['splits']['validation'][0]
  for r in records:
   if r['encoder']!=encoder:continue
   pred=run.predict(r,data[sample]);metrics=core.gene_metrics(data[sample]['y'][:,ix(r['panel'])],pred);stored=pd.read_csv(A/'logs/validation'/f'{r["run_id"]}.csv')
   for k,v in metrics.items():assert np.allclose(v,stored[k],rtol=1e-9,atol=1e-9,equal_nan=True),(r['run_id'],k)
   count+=1
 result={'timestamp':now(),'status':'PASS','models_replayed_on_validation':count,'test_targets_loaded':False,'selection_lock_sha256':sha(A/'SELECTION_LOCK.json')};dump(A/'logs/PRETEST_SELECTION_AUDIT.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):main()
