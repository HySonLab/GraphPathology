"""C1 numerical contracts, train-only scaling, metrics and compact CPU models."""
import copy
import os
import threading
import time
import numpy as np
from scipy.stats import rankdata
from sklearn.linear_model import Ridge
import torch
from torch import nn

class Scaler:
    def fit(self, values):
        values=np.asarray(values,dtype=np.float64)
        if values.ndim!=2 or not np.isfinite(values).all():
            raise ValueError('Expected a finite matrix')
        self.mean=values.mean(0); self.scale=values.std(0)
        self.scale[self.scale<1e-8]=1.
        return self
    def transform(self, values):
        return ((values-self.mean)/self.scale).astype(np.float32)
    def inverse(self, values):
        return np.asarray(values,dtype=np.float64)*self.scale+self.mean
    def save(self, path):
        np.savez_compressed(path,mean=self.mean,scale=self.scale)


def column_corr(a,b):
    a=np.asarray(a,dtype=np.float64);b=np.asarray(b,dtype=np.float64)
    if a.shape!=b.shape or a.ndim!=2 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Correlation matrices must have matching finite shapes')
    ac=a-a.mean(0);bc=b-b.mean(0)
    # Detect constant columns by range, not round-off in mean subtraction.
    valid=(np.ptp(a,axis=0)>1e-12)&(np.ptp(b,axis=0)>1e-12)
    out=np.full(a.shape[1],np.nan)
    denom=np.sqrt((ac*ac).sum(0)*(bc*bc).sum(0))
    valid &= denom>0
    out[valid]=np.clip((ac*bc).sum(0)[valid]/denom[valid],-1,1)
    return out


def gene_metrics(y,p):
    y=np.asarray(y,dtype=float);p=np.asarray(p,dtype=float)
    pearson=column_corr(y,p)
    spearman=column_corr(rankdata(y,axis=0),rankdata(p,axis=0))
    return {'pearson':pearson,'spearman':spearman,'mse':np.mean((y-p)**2,axis=0),
            'mae':np.mean(abs(y-p),axis=0),'target_variance':y.var(0),
            'prediction_variance':p.var(0)}


def finite_stat(x,method='mean'):
    x=np.asarray(x);x=x[np.isfinite(x)]
    return float(getattr(np,method)(x)) if len(x) else float('nan')


def summarize(metrics):
    out={}
    for key in ['pearson','spearman','mse','mae']:
        for stat in ['mean','median']:
            out[f'{stat}_{key}']=finite_stat(metrics[key],stat)
    p=metrics['pearson'];out['n_genes']=len(p);out['n_defined_pearson']=int(np.isfinite(p).sum())
    for threshold in [0,.1,.2,.3]:
        tag=str(threshold)
        out[f'n_pearson_gt_{tag}']=int(np.sum(p>threshold))
        out[f'fraction_pearson_gt_{tag}']=float(np.mean(p>threshold))
    out['median_prediction_variance']=float(np.median(metrics['prediction_variance']))
    return out


def correlation_matrix(x):
    x=np.asarray(x,dtype=float)
    centered=x-x.mean(0);norm=np.linalg.norm(centered,axis=0)
    valid=(np.ptp(x,axis=0)>1e-12)&(norm>0)
    z=np.zeros_like(centered);z[:,valid]=centered[:,valid]/norm[valid]
    matrix=np.clip(z.T@z,-1,1)
    matrix[~valid,:]=np.nan;matrix[:,~valid]=np.nan
    return matrix


def coexpression(y,p):
    truth,pred=correlation_matrix(y),correlation_matrix(p)
    index=np.triu_indices(y.shape[1],1)
    a,b=truth[index],pred[index];valid=np.isfinite(a)&np.isfinite(b)
    r=float(column_corr(a[valid,None],b[valid,None])[0]) if valid.sum()>1 else float('nan')
    error=float(np.sqrt(np.mean((a[valid]-b[valid])**2))) if valid.any() else float('nan')
    return {'coexpression_pearson':r,'coexpression_pair_rmse':error,
            'n_valid_pairs':int(valid.sum()),'n_total_pairs':len(a)},truth,pred

class PeakRSS:
    """Sample process RSS, including data held before the fit begins."""
    def __enter__(self):
        self.peak=0;self.stop=threading.Event()
        def sample():
            while not self.stop.is_set():
                try:
                    with open('/proc/self/statm') as f:
                        self.peak=max(self.peak,int(f.read().split()[1])*os.sysconf('SC_PAGE_SIZE'))
                except OSError: pass
                self.stop.wait(.02)
        self.thread=threading.Thread(target=sample,daemon=True);self.thread.start()
        return self
    def __exit__(self,*args):
        self.stop.set();self.thread.join()


def mlp_model(d,g,config):
    return nn.Sequential(nn.Linear(d,config['hidden']),nn.GELU(),nn.Dropout(config['dropout']),nn.Linear(config['hidden'],g))


def train_mlp(x,y,xv,yv,config,training,seed,fixed_epochs=None):
    torch.manual_seed(seed)
    model=mlp_model(x.shape[1],y.shape[1],config)
    optimizer=torch.optim.AdamW(model.parameters(),lr=training['learning_rate'],weight_decay=training['weight_decay'])
    xt,yt=torch.from_numpy(x),torch.from_numpy(y)
    xv_t=torch.from_numpy(xv) if xv is not None else None
    yv_t=torch.from_numpy(yv) if yv is not None else None
    history=[];best=float('inf');best_state=None;best_epoch=0;stale=0;train_s=0.;val_s=0.
    for epoch in range(1,(fixed_epochs or training['max_epochs'])+1):
        start=time.perf_counter();model.train();order=torch.randperm(len(x));loss_sum=0.
        for begin in range(0,len(x),training['batch_size']):
            idx=order[begin:begin+training['batch_size']]
            optimizer.zero_grad(set_to_none=True);loss=nn.functional.mse_loss(model(xt[idx]),yt[idx])
            if not torch.isfinite(loss): raise FloatingPointError('Nonfinite training loss')
            loss.backward();optimizer.step();loss_sum+=float(loss.detach())*len(idx)
        train_s+=time.perf_counter()-start
        if fixed_epochs is None:
            start=time.perf_counter();model.eval()
            with torch.inference_mode():
                total=0.
                for begin in range(0,len(xv_t),2048):
                    sl=slice(begin,begin+2048)
                    total+=float(nn.functional.mse_loss(model(xv_t[sl]),yv_t[sl],reduction='sum'))
                val=total/yv_t.numel()
            val_s+=time.perf_counter()-start
            if not np.isfinite(val): raise FloatingPointError('Nonfinite validation loss')
            if val<best-training['min_delta']:
                best=val;best_state=copy.deepcopy(model.state_dict());best_epoch=epoch;stale=0
            else: stale+=1
            history.append({'epoch':epoch,'training_mse':loss_sum/len(x),'validation_mse':val})
            if stale>=training['patience']: break
        else:
            history.append({'epoch':epoch,'training_mse':loss_sum/len(x)})
    if fixed_epochs is None: model.load_state_dict(best_state)
    else: best_epoch=fixed_epochs
    model.eval()
    return model,{'train_seconds':train_s,'validation_seconds':val_s,'best_epoch':best_epoch,
                  'epochs_run':epoch,'history':history,'parameters':sum(p.numel() for p in model.parameters())}


def predict(model,x,kind):
    if kind=='B0': return np.repeat(np.asarray(model)[None,:],len(x),axis=0)
    if kind=='B1':
        result=model.predict(x)
        return result[:,None] if result.ndim==1 else result
    with torch.inference_mode():
        return np.concatenate([model(torch.from_numpy(x[i:i+2048])).numpy() for i in range(0,len(x),2048)])


def permute_training(x,y,kind,seed):
    index=np.random.default_rng(seed).permutation(len(x))
    if kind=='NC1': return x,y[index],index
    if kind=='NC2': return x[index],y,index
    raise ValueError('Unknown negative control')
