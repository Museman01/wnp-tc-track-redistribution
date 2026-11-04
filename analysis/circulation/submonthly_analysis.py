"""Storm-matched steering skill, circulation statistics and path decomposition."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from core.paths import PROJECT
from analysis.circulation import submonthly_statistics as s


def wind_statistics(root, output, inference=False):
    base=root/'data/processed/submonthly_steering'
    frame=pd.read_csv(base/'matched_winds.csv',parse_dates=['iso_time'],float_precision='round_trip')
    reference=pd.read_csv(root/'results/circulation/submonthly_steering/wind_period_comparisons.csv',float_precision='round_trip')
    rows=[]
    keys=['geometry','layer','window_start','temporal_resolution','family']
    for group, expected in reference.groupby(keys,sort=False):
        geometry,layer,start,resolution,family=group
        part=frame[frame.season>=start]
        if geometry=='annulus_5_7':part=part[part.iso_time.dt.month.between(6,10)]
        years=np.arange(start,2026);half=len(years)//2
        prefix=geometry+'_'+layer
        if resolution=='sixhour_advantage_over_sameyear_month':
            high=s.annual_sufficient(part,[prefix+'_u',prefix+'_v'],years)
            month=s.annual_sufficient(part,[prefix+'_u_monthmatched',prefix+'_v_monthmatched'],years)
            def metric(ix):
                a=s.pooled_metric(*high,ix);b=s.pooled_metric(*month,ix)
                return np.array([a[2]-b[2],b[5]-a[5]])
            names=['vector_correlation_advantage','vector_RMSE_reduction']
        else:
            suffix='_monthmatched' if resolution=='same_year_month_synoptic_mean' else ''
            sufficient=s.annual_sufficient(part,[prefix+'_u'+suffix,prefix+'_v'+suffix],years)
            def metric(ix):return s.pooled_metric(*sufficient,ix)[:5]
            names=s.METRICS[:5]
        if inference:result=s.infer(metric,years,undefined_bootstrap='report_na')
        else:
            early=metric(np.arange(half));late=metric(np.arange(half,len(years)))
            result={'early':early,'late':late,'difference':late-early}
        for i,name in enumerate(names):
            early_n=int((part.season<start+half).sum());late_n=int((part.season>=start+half).sum())
            rows.append(dict(zip(keys,group),metric=name,n_segments_early=early_n,n_segments_late=late_n,
                             **{k:float(result[k][i]) for k in ['early','late','difference','ci_low','ci_high','p'] if k in result}))
    out=pd.DataFrame(rows)
    if inference:
        out['q']=np.nan
        for _,indices in out.groupby('family',sort=False).groups.items():
            values=out.loc[indices,'p'].to_numpy();valid=np.isfinite(values)
            adjusted=s.bh_fdr(np.where(valid,values,1.));adjusted[~valid]=np.nan
            out.loc[indices,'q']=adjusted
    out.to_csv(output/'wind_period_comparisons.csv',index=False)
    return out


def regime_statistics(root, output, inference=False):
    base=root/'data/processed/submonthly_steering'
    daily=pd.read_csv(base/'daily_regimes.csv',parse_dates=['date'],float_precision='round_trip')
    annual=s.regime_annual(daily,4)
    fixed=[]
    for year in s.YEARS:
        with np.load(base/f'daily_winds/{year}.npz') as z:values=z['fixed']
        fixed.append(dict(year=year,reference_u_mean=values[:,0].mean(),reference_v_mean=values[:,1].mean(),reference_eastward_fraction=values[:,2].mean()))
    annual=annual.merge(pd.DataFrame(fixed),on='year',validate='one_to_one')
    annual.to_csv(output/'annual_regimes.csv',index=False)
    reference=pd.read_csv(root/'results/circulation/submonthly_steering/regime_period_comparisons.csv')
    rows=[]
    for (start,family),group in reference.groupby(['window_start','family'],sort=False):
        part=annual[annual.year>=start];years=np.arange(start,2026);half=len(years)//2;subset=[]
        for name in group.metric:
            values=part[name].to_numpy()
            def metric(ix):return np.array([np.nanmean(values[ix])])
            if inference:result=s.infer(metric,years,undefined_bootstrap='report_na')
            else:
                early=metric(np.arange(half));late=metric(np.arange(half,len(years)))
                result={'early':early,'late':late,'difference':late-early}
            subset.append(dict(metric=name,window_start=start,window_end=2025,family=family,
                               **{k:float(result[k][0]) for k in ['early','late','difference','ci_low','ci_high','p'] if k in result}))
        if inference:
            values=np.array([row['p'] for row in subset]);valid=np.isfinite(values)
            adjusted=s.bh_fdr(np.where(valid,values,1.));adjusted[~valid]=np.nan
            for row,q in zip(subset,adjusted):row['q']=q
        rows.extend(subset)
    out=pd.DataFrame(rows);out.to_csv(output/'regime_period_comparisons.csv',index=False)
    return out


def joint_fields(root):
    base=root/'data/processed/submonthly_steering'
    tracks=pd.read_csv(base/'domain_ts_standard.csv',float_precision='round_trip')
    labels=pd.read_csv(base/'daily_regimes.csv',parse_dates=['date'])
    labels['date']=pd.to_datetime(labels.date,utc=True)
    points=s.label_points(tracks,labels.set_index('date').regime)
    joint,_=s.path_joint(points,list(range(4)))
    return joint


def path_statistics(root, output, inference=False):
    joint=joint_fields(root)
    if inference:result,summary,_,_,_=s.inference(joint)
    else:
        result=s.decompose(joint[:30].mean(0),joint[30:].mean(0))
        summary=pd.DataFrame([dict(component=k,projection_percent=s.project(result[k],result['delta'])) for k in ['between','within']])
    np.savez_compressed(output/'path_decomposition.npz',observed=result['delta'],frequency=result['between'],within=result['within'])
    summary.to_csv(output/'path_decomposition_inference.csv',index=False)
    return result,summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis',choices=['wind','regimes','path','all'])
    parser.add_argument('--root',type=Path,default=PROJECT)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--inference',action='store_true',help='Use 9999 block permutations and 4999 bootstrap replicates')
    args=parser.parse_args();root=args.root.resolve();output=args.output.resolve()
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    if args.analysis in ['wind','all']:wind_statistics(root,output,args.inference)
    if args.analysis in ['regimes','all']:regime_statistics(root,output,args.inference)
    if args.analysis in ['path','all']:path_statistics(root,output,args.inference)
    print(str(output))


if __name__=='__main__':main()
