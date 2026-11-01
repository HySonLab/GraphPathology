import numpy as np
import pytest
from graphpathology.extension_evaluate import score_set,cohort


def test_same_H_image_metrics_exact_and_all_unobserved_corruption():
 rng=np.random.default_rng(83);n=80;d={'targets':rng.normal(size=(n,50)),'xy_um':rng.normal(size=(n,2))*200,'genes':np.array([f'g{i}' for i in range(50)])};image=rng.normal(size=(n,50));hidden=np.arange(60,80);cache={};item={'id':'synthetic','origin':'new','variant':'matched_local','seed':0}
 base={'design':'B2_identical_H','condition':'distributed','budget':.2,'mask_seed':0,'geometry':1,'observed':np.arange(16),'hidden':hidden,'mask_sha256':'test'}
 a=score_set(d,image,base,item,'TENX99');b=score_set(d,image,base,item,'TENX99',cache);c=score_set(d,image,{**base,'observed':np.arange(20,36),'condition':'region_exclusion'},item,'TENX99',cache)
 assert a==b
 for metric in ['median_pearson','mean_mse','mean_mae']:assert a[0][0][metric]==c[0][0][metric]
 assert a[-1]==b[-1]==c[-1]==4
 assert len(cache)==1
 assert cohort('BASS_P01')=='Bassiouni' and cohort('CID4290')=='Wu'
