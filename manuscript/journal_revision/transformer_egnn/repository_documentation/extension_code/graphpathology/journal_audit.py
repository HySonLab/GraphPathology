"""All available original-reference numerical checks for measured Wu replay outputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .assets import sha


def audit_wu(root,output):
    source=root/'experiments/journal_revision/results/wu_anchors';old=root/'checkpoints/legacy/C9/C9B'
    panels=pd.concat([pd.read_csv(p) for p in source.glob('*/PANEL_METRICS.csv.gz')],ignore_index=True)
    selected=((panels.regime=='V1')&(panels.method=='idw_k16_p2'))|((panels.regime=='V2')&(panels.method=='uniform_k16'))
    actual=panels[selected];expected=pd.read_csv(old/'C9B_PANEL_METRICS.csv')
    expected=expected[(expected.method=='M1')&(expected.strategy=='random')&expected.control.fillna('').eq('')]
    keys=['patient','regime','geometry','seed','budget'];expected=expected.rename(columns={'patient_id':'patient'})
    metrics=['mean_pearson','median_pearson','mean_mse','median_mse','mean_mae','median_mae']
    merged=actual[keys+metrics].merge(expected[keys+metrics],on=keys,suffixes=('_actual','_expected'),validate='one_to_one',how='outer',indicator=True)
    if len(merged)!=360 or not merged['_merge'].eq('both').all():raise ValueError('Expected all 360 original random-reference runs')
    errors={}
    for metric in metrics:
        a=merged[metric+'_actual'];b=merged[metric+'_expected'];np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-10,equal_nan=True);errors[metric]=float(np.nanmax(abs(a-b)))
    actual=pd.concat([pd.read_csv(p) for p in source.glob('*/GENE_METRICS.csv.gz')],ignore_index=True)
    actual=actual[((actual.regime=='V1')&(actual.method=='idw_k16_p2'))|((actual.regime=='V2')&(actual.method=='uniform_k16'))]
    expected=pd.read_csv(old/'C9B_GENE_METRICS.csv');expected=expected[(expected.method=='M1')&(expected.strategy=='random')&expected.control.fillna('').eq('')].rename(columns={'patient_id':'patient'})
    keys=keys+['gene'];metrics=['pearson','mse','mae']
    m=actual[keys+metrics].merge(expected[keys+metrics],on=keys,suffixes=('_actual','_expected'),validate='one_to_one',how='outer',indicator=True)
    if len(m)!=18000 or not m['_merge'].eq('both').all():raise ValueError('Expected all 18,000 original-reference gene rows')
    gene_errors={}
    for metric in metrics:
        a=m[metric+'_actual'];b=m[metric+'_expected'];np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-10,equal_nan=True);gene_errors[metric]=float(np.nanmax(abs(a-b)))
    result={'status':'PASS_PRIMARY_ALL_ORIGINAL_RANDOM_REFERENCES','run_rows':len(merged),'gene_rows':len(m),'primary_gene_values_checked':len(m)*3,
            'summary_max_abs_error':errors,'gene_max_abs_error':gene_errors,'atol':1e-10,'rtol':1e-10,
            'secondary_spearman':'Not covered by exactness claim; flagged near-tied rank sensitivity',
            'original_panel_sha256':sha(old/'C9B_PANEL_METRICS.csv'),'original_gene_sha256':sha(old/'C9B_GENE_METRICS.csv')}
    if output.exists():raise FileExistsError(output)
    output.write_text(json.dumps(result,indent=2)+'\n');return result
