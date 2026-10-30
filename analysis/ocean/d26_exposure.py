"""P5: bounded primary-record D26 decomposition with a matched OHC300 comparator."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import sys

import numpy as np

import pandas as pd

import xarray as xr

import argparse

from core.paths import PROJECT,output_path,write_json,import_source



def main():

    present=[y for y in range(1966,2026) if output_path(f'data/derived/d26/oras5_d26_{y}.nc').is_file()]

    if len(present)!=60:raise RuntimeError(f'D26 incomplete: {len(present)}/60 years; no partial period comparison allowed')

    src=_LAYOUT_ROOT / 'analysis/ocean'

    sys.path.insert(0,str(src));dec=import_source('ohc_pure_decomposition','analysis/ocean/ohc_decomposition.py')

    matched=pd.read_csv(_LAYOUT_ROOT / 'results/ocean/ohc300/matched_points.csv.gz')

    matched=matched.loc[matched.agency.eq('PRIMARY')&matched.definition.eq('ts_only')&matched.match_status.eq('matched')].copy()

    keys=dec.state_table(matched)

    full=dec.attach_key_ids(matched,keys)

    months=keys.ohc_month.to_numpy(int)-1;gy=keys.grid_y.to_numpy(int);gx=keys.grid_x.to_numpy(int)

    D=np.empty((60,len(keys)));H=np.empty_like(D);geometry=True

    for row,year in enumerate(range(1966,2026)):

        with xr.open_dataset(output_path(f'data/derived/d26/oras5_d26_{year}.nc')) as d, xr.open_dataset(PROJECT/f'data/derived/ohc300/oras5_ohc300_{year}.nc') as h:

            if not (d.d26.attrs['units']=='m'):raise ValueError('D26 units')

            same=d.nav_lat.shape==h.nav_lat.shape and np.allclose(d.nav_lat,h.nav_lat,atol=5e-5,rtol=0) and np.allclose(d.nav_lon%360,h.nav_lon%360,atol=5e-5,rtol=0)

            if not same:raise ValueError('D26/OHC coordinate mismatch; spatial rematching required')

            D[row]=d.d26.values[months,gy,gx];H[row]=h.ohc300.values[months,gy,gx]

    support=np.isfinite(D).all(axis=0)&np.isfinite(H).all(axis=0)

    bad=set(keys.loc[~support,'key_id'].astype(int))

    common=full.loc[~full.key_id.isin(bad)].copy()

    # Reindex keys after any complete-support exclusion, preserving the matching rule.

    keymap={old:new for new,old in enumerate(keys.loc[support,'key_id'])}

    common['key_id']=common.key_id.map(keymap).astype(int)

    D=D[:,support];H=H[:,support]

    keptkeys=keys.loc[support].copy();keptkeys['key_id']=np.arange(len(keptkeys))

    keptkeys.to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_state_keys.csv',index=False)

    coverage=[]

    for year in range(1966,2026):

        a=full.loc[full.ohc_year.eq(year)];b=common.loc[common.ohc_year.eq(year)]

        coverage.append({'year':year,'original_matched_points':len(a),'common_support_points':len(b),

            'common_support_point_fraction':len(b)/len(a),'native_zero_state_fraction':float(np.mean(D[year-1966]==0))})

    pd.DataFrame(coverage).to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_coverage.csv',index=False)

    summaries=[];components=[];bootstraps=[];annual=[];checks=[]

    for land in ['ocean_points_only','nearest_ocean_all_points']:

        view=dec.analysis_view(common,'calendar_year',land,1966,2025)

        W,storms,points,effective=dec.annual_weight_matrix(view,list(range(1966,2026)),len(keptkeys),'storm_normalized_equal_year')

        for variable,F,units in [('D26',D,'m'),('OHC300_same_support',H,'J m-2')]:

            for comparison,p1,p2 in [('1966_1995_vs_1996_2025',(1966,1995),(1996,2025)),('1982_2003_vs_2004_2025',(1982,2003),(2004,2025))]:

                i=np.arange(p1[0]-1966,p1[1]-1966+1);j=np.arange(p2[0]-1966,p2[1]-1966+1)

                result,arrays=dec.symmetric_decomposition(W[i],F[i],W[j],F[j],1e-8 if units=='m' else 1e-4,1e-12)

                ident={'variable':variable,'units':units,'land_treatment':land,'comparison':comparison,

                    'agency':'PRIMARY','definition':'ts_only','weighting':'storm_normalized_equal_year','year_assignment':'calendar_year'}

                summaries.append({**ident,**result})

                bs=dec.bootstrap_summary(W[i],F[i],W[j],F[j],3,1999,20260801,1e-8 if units=='m' else 1e-4,1e-12)

                for key,value in ident.items():bs[key]=value

                bootstraps.append(bs)

                if variable=='D26':

                    cell=keptkeys[['key_id','ohc_month','grid_lat','grid_lon']].copy()

                    for key in ['redistribution_cell','ocean_cell','weight_p1','weight_p2']:cell[key]=arrays[key]

                    for key,value in ident.items():cell[key]=value

                    components.append(cell)

                checks.append(abs(result['observed_closure_residual']))

                print('Parallel exposure',variable,land,comparison,'delta',result['delta_exposure'],flush=True)

            exp=np.einsum('ij,ij->i',W,F)

            fixed=F@W.mean(axis=0)

            for y,value,n,background in zip(range(1966,2026),exp,points,fixed):annual.append({**ident,'comparison':'all_years','year':y,'annual_exposure':value,'fixed_path_background_exposure':background,'n_points':int(n)})

    summary=pd.DataFrame(summaries);summary.to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_parallel_summary.csv',index=False)

    pd.concat(bootstraps,ignore_index=True).to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_parallel_bootstrap.csv',index=False)

    pd.concat(components,ignore_index=True).to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_cell_components.csv.gz',index=False,compression='gzip')

    ann=pd.DataFrame(annual);ann.to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_annual_exposure.csv',index=False)

    # Reproduce original OHC point estimates whenever no support was removed.

    reference=pd.read_csv(_LAYOUT_ROOT / 'results/ocean/ohc300/decomposition_summary.csv')

    errors=[]

    for row in summary.loc[summary.variable.eq('OHC300_same_support')].itertuples():

        r=reference.loc[reference.agency.eq('PRIMARY')&reference.definition.eq('ts_only')&reference.year_assignment.eq('calendar_year')&reference.land_treatment.eq(row.land_treatment)&reference.weighting.eq('storm_normalized_equal_year')&reference.comparison.eq(row.comparison)]

        for name in ['delta_exposure','redistribution_component','ocean_component','covariability_component']:

            old=r[name]

            if len(old)!=1:raise ValueError('Original OHC reference not unique')

            errors.append(float(getattr(row,name)-old.iloc[0]))

    maxerr=max(map(abs,errors))

    if not bad and maxerr>1e-3:raise AssertionError(f'OHC comparator does not reproduce original: {maxerr}')

    # Product seam is descriptive; no numerical adjustment or causal claim.

    seam=ann.loc[ann.variable.eq('D26')].pivot(index='year',columns='land_treatment',values='annual_exposure')

    seam.to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_product_seam_annual.csv')

    seamrows=[]

    for land,g in ann.loc[ann.variable.eq('D26')].groupby('land_treatment'):

        g=g.set_index('year').sort_index();b=g.fixed_path_background_exposure;diff=b.diff().dropna()

        seamrows.append({'land_treatment':land,'fixed_path_2015_minus_2014_m':float(diff.loc[2015]),

            'absolute_jump_rank_of_59':int((abs(diff)>=abs(diff.loc[2015])).sum()),

            'background_2015_2019_minus_2010_2014_m':float(b.loc[2015:2019].mean()-b.loc[2010:2014].mean()),

            'diagnostic':'descriptive fixed pooled path-month weights; production change is confounded with climate variability; no correction'})

    pd.DataFrame(seamrows).to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_product_seam_diagnostic.csv',index=False)

    write_json('results/diagnostics/d26_parallel.json',{'status':'PASS','configurations':'one record × one TS definition × one year assignment × one weighting × two land treatments × two time windows; no 64-way matrix',

        'n_state_keys_before':len(keys),'n_state_keys_retained':len(keptkeys),'removed_state_keys':len(bad),

        'coordinate_equality':True,'zero_handling':'native delivered zero, never an imputation for land or missing data',

        'maximum_absolute_closure_residual_mixed_units':max(checks),'ohc_reference_maxabs_j_m2':maxerr,

        'uncertainty':'1999 circular-three-year block bootstrap draws, paired annual occupancy and field within each period',

        'production_seam':'consolidated through 2014, operational from 2015; retained without correction',

        'limits':['Monthly background exposure only, not prestorm structure or storm-induced feedback.','1966-1978 is a backward spin-up extension.',

                  '2014/2015 seam is not experimentally isolated from natural variability.']})



def seam_only():

    keys=pd.read_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_state_keys.csv').sort_values('key_id')

    components=pd.read_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_cell_components.csv.gz')

    components=components.loc[components.comparison.eq('1966_1995_vs_1996_2025')]

    pools={land:.5*(g.set_index('key_id').loc[keys.key_id,'weight_p1'].to_numpy()+g.set_index('key_id').loc[keys.key_id,'weight_p2'].to_numpy()) for land,g in components.groupby('land_treatment')}

    values={land:[] for land in pools}

    for year in range(1966,2026):

        with xr.open_dataset(output_path(f'data/derived/d26/oras5_d26_{year}.nc')) as ds:

            f=ds.d26.values[keys.ohc_month.to_numpy(int)-1,keys.grid_y.to_numpy(int),keys.grid_x.to_numpy(int)]

        for land,pool in pools.items():values[land].append(float(pool@f))

    frame=pd.DataFrame(values,index=range(1966,2026));frame.index.name='year'

    frame.to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_fixed_path_background_annual.csv')

    rows=[]

    for land in frame:

        b=frame[land];diff=b.diff().dropna()

        rows.append({'land_treatment':land,'fixed_path_2015_minus_2014_m':float(diff.loc[2015]),

            'absolute_jump_rank_of_59':int((abs(diff)>=abs(diff.loc[2015])).sum()),

            'background_2015_2019_minus_2010_2014_m':float(b.loc[2015:2019].mean()-b.loc[2010:2014].mean()),

            'diagnostic':'descriptive pooled fixed path-month weights; production change confounded with climate variability; no correction'})

    pd.DataFrame(rows).to_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_product_seam_diagnostic.csv',index=False)

    print('Fixed-path D26 production-seam diagnostics saved.')



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    parser=argparse.ArgumentParser();parser.add_argument('--seam-only',action='store_true');args=parser.parse_args()

    if args.seam_only:seam_only()

    else:main()
