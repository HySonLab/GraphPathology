import io
import pandas as pd
import numpy as np
from graphpathology import csvio as runner
def test_single_control_csv_preserves_identity_and_values():
 blocks=[]
 for control in ['TRUE','MATCHED','MEAN','LOCAL']:
  source=f'control,gene,pearson,mse\n{control},GENE1,0.123456789,1.23456789\n'
  d=runner.typed_read_csv(io.StringIO(source));assert d.control.iloc[0]==control
  np.testing.assert_array_equal(d[['pearson','mse']].to_numpy(),[[.123456789,1.23456789]])
  blocks.append(d)
 pivot=pd.concat(blocks).pivot(index='gene',columns='control',values='pearson')
 assert set(pivot.columns)=={'TRUE','MATCHED','MEAN','LOCAL'}
 assert float((pivot['TRUE']-pivot['MATCHED']).iloc[0])==0
