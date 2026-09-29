"""Run the immutable C4B analysis with explicit text typing for CSV control labels.

Pandas infers a single-label TRUE column as Boolean in individual gene-metric
files. This adapter restores the JSON-recorded text label; scientific formulas
and the locked analysis source remain unchanged.
"""
from common import *
original=pd.read_csv
def typed_read_csv(*args,**kwargs):
 dtype=kwargs.pop('dtype',{}) or {};kwargs['dtype']={**dtype,'control':str};return original(*args,**kwargs)
if __name__=='__main__':
 pd.read_csv=typed_read_csv
 spec=importlib.util.spec_from_file_location('c4b_locked_analysis',OUT/'b_analyze.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 try:module.main()
 finally:pd.read_csv=original
