"""Generate Tables S13/S14 and the input arrays for Figures S20/S21."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from core.paths import PROJECT
from analysis.circulation.submonthly_analysis import joint_fields
from analysis.circulation import submonthly_statistics as s


def products(root, output, refit=False):
    results=root/'results/circulation/submonthly_steering'
    winds=pd.read_csv(results/'wind_period_comparisons.csv',float_precision='round_trip')
    columns=['geometry','temporal_resolution','metric','n_segments_early','n_segments_late','early','late','difference','ci_low','ci_high','p','q']
    selected=winds[(winds.window_start==1966)&(winds.layer=='deep')&winds.metric.isin(s.METRICS[:5])]
    selected[columns].to_csv(output/'table_s13_motion_environment.csv',index=False)
    regime=pd.read_csv(results/'regime_period_comparisons.csv',float_precision='round_trip')
    regime[['metric','window_start','window_end','early','late','difference','ci_low','ci_high','p','q','family']].to_csv(output/'table_s14_regimes.csv',index=False)
    base=root/'data/processed/submonthly_steering'
    with np.load(base/'circulation_model.npz') as model:latitude=model['latitude'];longitude=model['longitude']
    fields=np.stack([np.load(base/f'daily_winds/{year}.npz')['field'] for year in s.YEARS])
    anomalies=(fields-fields.mean(0)).reshape(9180,2,18,32)
    labels=pd.read_csv(base/'daily_regimes.csv',parse_dates=['date'])
    composites=np.stack([anomalies[labels.regime.values==k].mean(0) for k in range(4)])
    annual=pd.read_csv(base/'annual_regimes.csv',float_precision='round_trip')
    primary=regime[regime.window_start==1966].set_index('metric')
    np.savez(output/'figS20_regimes.npz',latitude=latitude,longitude=longitude,composite_uv_anomaly=composites,years=annual.year.values,
             annual_frequency=annual[[f'frequency_{k}' for k in range(4)]].values,
             early_frequency=[primary.loc[f'frequency_{k}','early'] for k in range(4)],late_frequency=[primary.loc[f'frequency_{k}','late'] for k in range(4)])
    joint=joint_fields(root);result=s.decompose(joint[:30].mean(0),joint[30:].mean(0))
    np.savez(output/'figS21_decomposition.npz',lat_edges=np.arange(0,40.01,2.5),lon_edges=np.arange(100,180.01,2.5),
             observed_pp=result['delta']*100,frequency_pp=result['between']*100,within_pp=result['within']*100)
    if refit:
        from sklearn.decomposition import PCA
        from sklearn.cluster import KMeans
        from threadpoolctl import threadpool_limits
        points=pd.read_csv(base/'domain_ts_standard.csv',parse_dates=['iso_time'],float_precision='round_trip')
        points=points[points.iso_time.dt.month.between(6,10)].copy();points['date']=points.iso_time.dt.normalize()
        matrix=(anomalies*np.sqrt(np.cos(np.deg2rad(latitude)))[None,None,:,None]).reshape(9180,-1)
        rows=[dict(k=4,labels='published classification',frequency_percent=s.project(result['between'],result['delta']),within_percent=s.project(result['within'],result['delta']))]
        with threadpool_limits(limits=1):
            pcs=PCA(17,svd_solver='full').fit_transform(matrix)
            for k in range(4,9):
                assignment=KMeans(k,n_init=30,random_state=1).fit_predict(pcs)
                joined=points.merge(labels[['date']].assign(label=assignment),on='date',how='left',validate='many_to_one')
                values=np.zeros((60,k,16,32))
                for i,year in enumerate(s.YEARS):
                    group=joined[joined.year==year]
                    for index in range(k):
                        subset=group[group.label==index];values[i,index]=np.histogram2d(subset.lat,subset.lon,bins=[np.arange(0,40.01,2.5),np.arange(100,180.01,2.5)])[0]
                    values[i]/=values[i].sum()
                partition=s.decompose(values[:30].mean(0),values[30:].mean(0))
                rows.append(dict(k=k,labels='refit k-means (n_init=30, random_state=1)',frequency_percent=s.project(partition['between'],partition['delta']),within_percent=s.project(partition['within'],partition['delta'])))
        pd.DataFrame(rows).to_csv(output/'cluster_count_sensitivity.csv',index=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,default=PROJECT);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--refit-clusters',action='store_true',help='Refit k=4 through 8 with 17 principal components and seed 1')
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    args.output.mkdir(parents=True);products(args.root.resolve(),args.output,args.refit_clusters);print(str(args.output))


if __name__=='__main__':main()
