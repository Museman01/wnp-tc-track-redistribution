"""Storm anchors and period effects in recurvature, with circulation-regime adjustment.

The period-only and period-plus-regime models use season-clustered covariance,
4,999 within-period three-year block resamples and 9,999 three-year block
permutations. Seed 202406. Coefficients describe associations, not causality.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize
from core.paths import PROJECT
from analysis.tracks.pathway_closure import block_indices
from analysis.tracks.track_metrics import recurvature,tropical_stage
from analysis.shared.block_resampling import block_order,bh_fdr
YEARS=np.arange(1966,2026);EARLY=np.arange(30);LATE=np.arange(30,60);SEED=202406

def validate_labels(frame):
    f=frame.copy();f['date']=pd.to_datetime(f.date,utc=True).dt.normalize()
    dates=pd.date_range('1966-01-01','2025-12-31',tz='UTC');dates=dates[dates.month.isin(range(6,11))]
    if f.date.duplicated().any() or set(f.date)!=set(dates):raise ValueError('require exactly9180 seasonal dates')
    if not np.isfinite(f.regime).all() or not np.equal(f.regime,np.floor(f.regime)).all():raise ValueError('invalid regime label')
    keys=sorted(f.regime.unique().astype(int))
    if not 4<=len(keys)<=6:raise ValueError('expected4–6 regimes')
    return f.set_index('date').regime.astype(int),keys

def permutation_order(n,rng):
    if n%3:raise ValueError('year count must be divisible by3')
    return block_order(n,3,rng)

def storm_anchors(full,eligible,labels):
    out=[];excluded=[];full=full.copy();full['iso_time']=pd.to_datetime(full.iso_time,utc=True)
    for sid,g in full[full.sid.isin(eligible)].groupby('sid'):
        g=tropical_stage(g,34).reset_index(drop=True)
        if g.empty:raise ValueError('missing full TS life')
        north=np.flatnonzero(g.lat.to_numpy()>=20)
        ix=int(north[0]) if len(north) else int(np.argmin(np.unwrap(np.deg2rad(g.lon.to_numpy()))))
        a=g.iloc[ix];rec,rt,_,_=recurvature(g,34)
        row=dict(sid=sid,season=int(a.season),anchor_time=a.iso_time,anchor_method='first_report_ge20' if len(north) else 'earliest_unwrapped_westmost',recurving=int(rec),recurve_time=rt,full_ts_reports=len(g))
        if not 6<=a.iso_time.month<=10:excluded.append(dict(**row,reason='anchor_outside_June_October'));continue
        date=a.iso_time.normalize()
        if date not in labels:raise ValueError('missing anchor regime')
        out.append(dict(**row,regime=int(labels.loc[date])))
    if set(eligible)-set(full.sid):raise ValueError('eligible SID absent from full life catalogue')
    return pd.DataFrame(out),pd.DataFrame(excluded)

def binomial_fit(x,y,n,clusters):
    if np.linalg.matrix_rank(x)<x.shape[1]:raise ValueError('logistic rank deficient')
    fun=lambda b:np.sum(n*np.logaddexp(0,x@b)-y*(x@b));jac=lambda b:x.T@(n*expit(x@b)-y)
    fit=minimize(fun,np.zeros(x.shape[1]),jac=jac,method='BFGS',options={'gtol':1e-7,'maxiter':300})
    if np.max(abs(jac(fit.x)))>1e-5 or np.max(abs(fit.x))>40:raise ValueError('logistic separation/nonconvergence')
    p=expit(x@fit.x);h=x.T@((n*p*(1-p))[:,None]*x)
    if np.linalg.cond(h)>1e12:raise ValueError('logistic singular information')
    bread=np.linalg.inv(h);score=x*(y-n*p)[:,None];groups=np.unique(clusters);scores=np.array([score[clusters==g].sum(0) for g in groups]);G=len(groups);N=n.sum();K=x.shape[1]
    if G<=1 or N<=K:raise ValueError('insufficient logistic clusters')
    cov=bread@(scores.T@scores)@bread*G/(G-1)*(N-1)/(N-K)
    return fit.x,np.sqrt(np.maximum(np.diag(cov),0))

def logistic_counts(a,keys):
    n=np.zeros((60,len(keys)));s=n.copy()
    for i,y in enumerate(YEARS):
        for k,r in enumerate(keys):
            g=a[(a.season==y)&(a.regime==r)];n[i,k]=len(g);s[i,k]=g.recurving.sum()
    return n,s

def logistic_pair(n,s,order=None):
    order=np.arange(60) if order is None else order;n=n[order];s=s[order];K=n.shape[1]
    period=np.repeat(np.r_[np.zeros(30),np.ones(30)],K);reg=np.tile(np.eye(K)[:,1:],(60,1));x=np.column_stack([np.ones(60*K),period,reg]);keep=n.ravel()>0;cl=np.repeat(np.arange(60),K)[keep]
    b0,se0=binomial_fit(x[keep,:2],s.ravel()[keep],n.ravel()[keep],cl);b1,se1=binomial_fit(x[keep],s.ravel()[keep],n.ravel()[keep],cl)
    return np.array([b0[1],b1[1],b1[1]-b0[1]]),[se0[1],se1[1],np.nan]

def logistic_inference(a,keys,nboot=4999,nperm=9999):
    n,s=logistic_counts(a,keys);beta,se=logistic_pair(n,s);rng=np.random.default_rng(SEED);boot=np.empty((nboot,3))
    for i in range(nboot):boot[i]=logistic_pair(n,s,np.r_[block_indices(EARLY,3,rng),block_indices(LATE,3,rng)])[0]
    rng=np.random.default_rng(SEED);exceed=np.zeros(3,int)
    for i in range(nperm):exceed+=abs(logistic_pair(n,s,permutation_order(60,rng))[0])>=abs(beta)
    ci=np.quantile(boot,[.025,.975],axis=0)
    return pd.DataFrame(dict(model=['period_only','period_plus_regime','adjusted_minus_unadjusted'],beta=beta,odds_ratio=np.exp(beta),ci_low=ci[0],ci_high=ci[1],or_ci_low=np.exp(ci[0]),or_ci_high=np.exp(ci[1]),year_cluster_se=se,permutation_p=(exceed+1)/(nperm+1))),boot


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=PROJECT);p.add_argument('--output',type=Path,required=True);p.add_argument('--make-anchors',action='store_true');p.add_argument('--inference',action='store_true');a=p.parse_args()
    root=a.root.resolve();a.output.mkdir(parents=True,exist_ok=False)
    base=root/'data/processed/submonthly_steering';results=root/'results/circulation/submonthly_steering'
    if a.make_anchors:
        labels,keys=validate_labels(pd.read_csv(base/'daily_regimes.csv'))
        domain=pd.read_csv(base/'domain_ts_standard.csv');full=pd.read_csv(base/'full_ts_standard.csv')
        anchors,excluded=storm_anchors(full,set(domain.sid),labels)
        anchors.to_csv(a.output/'storm_anchors_and_recurvature.csv',index=False);excluded.to_csv(a.output/'storm_anchor_exclusions.csv',index=False)
    else:anchors=pd.read_csv(results/'storm_anchors_and_recurvature.csv');keys=sorted(anchors.regime.unique())
    if a.inference:
        frame,_=logistic_inference(anchors,keys)
        family=np.r_[pd.read_csv(results/'shapley_projection_summary.csv').permutation_p,pd.read_csv(results/'overall_period_change_tests.csv').permutation_p,frame.permutation_p]
        if len(family)!=11:raise ValueError('Expected the four decomposition, four overall and three logistic tests')
        frame['q_bh_scalar_family']=bh_fdr(family)[-3:]
    else:
        n,s=logistic_counts(anchors,keys);beta,se=logistic_pair(n,s)
        frame=pd.DataFrame(dict(model=['period_only','period_plus_regime','adjusted_minus_unadjusted'],beta=beta,odds_ratio=np.exp(beta),year_cluster_se=se))
    frame.to_csv(a.output/'logistic_period_summary.csv',index=False)

if __name__=='__main__':main()
