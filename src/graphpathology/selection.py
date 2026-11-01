import numpy as np
import pandas as pd
def select(rows,key):
 df=pd.DataFrame(rows);cols=['median_pearson','common_mse'];patient=df.groupby([key,'sample_id','regime','budget'],as_index=False)[cols].mean().groupby([key,'sample_id'],as_index=False)[cols].mean();summary=patient.groupby(key)[cols].mean();best=summary.median_pearson.max();eligible=summary[summary.median_pearson>=best-.005];chosen=eligible.sort_values(['common_mse'],kind='stable').index[0];return chosen,summary.reset_index().to_dict('records')


def require_training_patients(requested, forbidden):
    if set(requested)&set(forbidden):raise ValueError('Held-out patient access forbidden before target read')
    return tuple(requested)
