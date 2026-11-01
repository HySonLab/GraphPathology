"""Independent numerical checks on generated views; no model inference or rendering."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd

def validate(root,out):
    checks=[]
    def check(name):checks.append(name)
    def close(actual,expected,name):
        np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-12,equal_nan=True);check(name)
    rows=json.loads((out/'TARGET_STATUS.json').read_text())['targets']
    assert len(rows)==29 and len({x['target'] for x in rows})==29
    assert sum(x['status']=='export_view_generated' for x in rows)==26;check('29 unique targets,26 explicitly partial export views,3 blocked')
    legacy=root/'checkpoints/legacy'
    # Copied immutable tables must retain their exact bytes, including compressed full per-gene data.
    for folder in ['tables','supp_tables']:
        for p in (legacy/'C10/MANUSCRIPT'/folder).glob('*.csv*'):
            assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256((out/folder/p.name).read_bytes()).digest();check('preserved '+p.name)
    # Current Table2 printed cells: independent baseline transcribed from supplied PDF Table2.
    baseline={
      'Development':[[.695,.9110,.7606],[.792,.2684,.3904],[.822,.2787,.4173]],
      'HER2ST zero-shot':[[.143,4.3665,1.9376],[.194,.0680,.1730],[.255,.3481,.5153]],
      'Wu fresh validation':[[.059,2.9155,1.5578],[.359,.0987,.2179],[.335,.1002,.2321]],
      'Bassiouni second replication':[[.075,1.9673,1.2120],[.416,.0923,.2033],[.407,.1022,.2288]]}
    d=pd.read_csv(out/'tables/Table2.csv')
    for cohort,expected in baseline.items():
        actual=d[d.cohort==cohort].set_index('method').loc[['M0','M1','M3'],['median_pearson','mean_mse','mean_mae']].to_numpy()
        close(np.round(actual[:,0],3),np.array(expected)[:,0],cohort+' printed Pearson')
        close(np.round(actual[:,1:],4),np.array(expected)[:,1:],cohort+' printed MSE/MAE')
    # Main external curves and paired effects against independent original patient summary CSVs.
    for target,cp in [('Figure4','C9'),('Figure5','C10')]:
        trace=json.loads((out/(target+'_plotted_values.json')).read_text())
        d=pd.read_csv(legacy/f'{cp}/{cp}B/{cp}B_PATIENT_METRICS.csv');d=d[(d.strategy=='random')&d.method.isin(['M0','M1','M3'])]
        for ax,metric in [(1,'median_pearson'),(2,'mean_mse')]:
            for line,method in zip(trace[ax]['lines'],['M0','M1','M3']):
                g=d[(d.method==method)&(d.regime=='V1')].groupby('budget')[metric].median()
                close(line['x'],g.index.to_numpy()*100,target+' '+method+' coverage')
                close(line['y'],g.to_numpy(),target+' '+method+' '+metric)
        p=d[(d.budget==.2)&(d.regime=='V1')].pivot(index='patient_id',columns='method',values='median_pearson')
        for collection,a,b in zip(trace[3]['collections'],['M1','M3'],['M0','M1']):
            close(np.sort(np.array(collection['offsets'])[:,1]),np.sort(p[a]-p[b]),target+' paired '+a+'-'+b)
    t8=pd.read_csv(out/'tables/TableS8.csv');assert len(t8)==9;check('TableS8 editable architecture schema and nine rows')
    # Original SFig6 has 4448 rows; retain undisplayed data and missing-value pattern.
    original=pd.read_csv(legacy/'C6/C6A/C6A_PANEL_METRICS.csv');generated=pd.read_csv(out/'data/SFig6_source.csv')
    assert len(original)==len(generated)==4448
    pd.testing.assert_frame_equal(original,generated);check('4448-row SFig6 retained,including nulls')
    for row in rows:
        if row['target'].startswith('Figure') and row['status']=='export_view_generated':
            trace=json.loads((out/(row['target']+'_plotted_values.json')).read_text());assert trace
            check(row['target']+' plotted value trace exists')
    result={'checks_passed':len(checks),'checks':checks,'rtol':1e-12,'atol':1e-12,'printed_cells':'Exact after manuscript decimal rounding','revised_baseline_verified':False,'limitations':'Revised full arrays/templates/source renderers absent; tests prove preserved exports and selected plotted/table values, not exact revised-display equivalence.'}
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    print(json.dumps(validate(args.root.resolve(),args.output.resolve()),indent=2))
