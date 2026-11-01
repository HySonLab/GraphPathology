import pandas as pd
def typed_read_csv(*args,**kwargs):
 dtype=kwargs.pop('dtype',{}) or {};kwargs['dtype']={**dtype,'control':str};return pd.read_csv(*args,**kwargs)
