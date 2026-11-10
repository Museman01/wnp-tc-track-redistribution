"""Propagate layer-specific corridor indices through the existing HAC model definitions."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

from pathlib import Path
import numpy as np

import pandas as pd

import statsmodels.api as sm

from scipy.signal import detrend

from core.paths import PROJECT,output_path

from analysis.shared.inference_stats import bh_fdr



def main(output=None):
    target=(PROJECT / (output or 'results/circulation/steering_hac_layer_models.csv')).resolve()
    if not target.is_relative_to(PROJECT.resolve()):raise ValueError('Output must remain within the project')
    if target.exists():raise FileExistsError(target)
    target.parent.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(_LAYOUT_ROOT / 'results/circulation/steering_corridor_annual.csv')

    oni=pd.read_csv(PROJECT/'data/raw/indices/oni.csv').rename(columns={'season':'year'})

    pdo=pd.read_csv(PROJECT/'data/raw/indices/pdo.csv').rename(columns={'season':'year'})

    a=np.load(output_path('data/derived/circulation_monthly_coarse.npz'))

    z=a['z500'];z=z-z.mean(2,keepdims=True);lat=a['latitude'];lon=a['longitude']

    la=(lat>=15)&(lat<=35);lo=(lon>=110)&(lon<=160)

    zz=[float(z[(a['year']==y)&np.isin(a['month'],[6,7,8,9,10])][:,la][:,:,lo].mean()) for y in range(1966,2026)]

    rows=[]

    models={'oni':['jas_oni'],'oni_pdo':['jas_oni','pdo'],'eddy_wnpsh_oni_pdo':['eddy_height','jas_oni','pdo'],

            'steering_oni_pdo':['u_ms','v_ms','jas_oni','pdo']}

    for layer,g in d.groupby('layer_hpa'):

        g=g.sort_values('year').merge(oni,on='year').merge(pdo,on='year');g['eddy_height']=zz

        for scale in ['raw','detrended']:

            for model,predictors in models.items():

                data=g[['index_oos_z',*predictors]].dropna().copy()

                if scale=='detrended':data[:]=detrend(data.to_numpy(),axis=0)

                data=(data-data.mean())/data.std(ddof=1)

                fit=sm.OLS(data.index_oos_z,sm.add_constant(data[predictors])).fit(cov_type='HAC',cov_kwds={'maxlags':3})

                for term in fit.params.index:

                    rows.append({'layer':layer,'timescale':scale,'model':model,'term':term,

                        'coefficient_standardized':float(fit.params[term]),'p_hac':float(fit.pvalues[term]),

                        'ci_low':float(fit.conf_int().loc[term,0]),'ci_high':float(fit.conf_int().loc[term,1]),'n':int(fit.nobs)})

    result=pd.DataFrame(rows);result['q_bh_within_layer_timescale']=np.nan

    for _,idx in result.loc[~result.term.eq('const')].groupby(['layer','timescale']).groups.items():

        result.loc[idx,'q_bh_within_layer_timescale']=bh_fdr(result.loc[idx,'p_hac'])

    result.to_csv(target,index=False)

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--execute',action='store_true');parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if not args.execute:parser.error('add --execute to calculate steering regressions')
    main(args.output)

