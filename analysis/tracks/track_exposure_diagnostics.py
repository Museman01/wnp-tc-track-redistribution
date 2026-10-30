'Reference track contrasts, model scenarios and absolute activity metrics.'

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import json

import numpy as np

import pandas as pd

from core.paths import PROJECT,output_path,write_json

from analysis.shared.agency_catalogues import all_catalogues, primary_ts, domain, annual_fields

from analysis.tracks.track_metrics import recurvature



def bam_report(primary):

    a=np.load(_LAYOUT_ROOT / 'results/models/counterfactual/fields.npz')

    fields=annual_fields(primary,np.arange(1966,2026))

    d=fields[30:].mean(0)-fields[:30].mean(0)

    archived=a['d_obs'].ravel()

    error=float(np.max(abs(d-archived)))

    if error>1e-12:raise AssertionError(f'Observed difference disagrees with archive: {error}')

    tv=lambda x:float(abs(x).sum()/2)

    scenarios=[]

    for name in ['d_full','d_steer','d_loc','d_life']:

        x=a[name].ravel()

        scenarios.append({'scenario':name,'spatial_r':float(np.corrcoef(x,d)[0,1]),'tv':tv(x),

            'tv_amplitude_ratio':tv(x)/tv(d),'projection_on_observed_axis':float(x@d/(d@d))})

    pd.DataFrame(scenarios).to_csv(_LAYOUT_ROOT / 'results/models/bam_existing_scenarios.csv',index=False)

    baseline=[]

    xe=np.arange(101.25,180,2.5); ye=np.arange(1.25,40,2.5)

    xx,yy=np.meshgrid(xe,ye)

    for label,key,obs in [('early','f_BASE_P1',fields[:30].mean(0)),('late','f_BASE_P2',fields[30:].mean(0))]:

        sim=a[key].ravel()

        baseline.append({'period':label,'r':float(np.corrcoef(obs,sim)[0,1]),'tv_error':tv(sim-obs),

            'observed_share_north25':float(obs[yy.ravel()>=25].sum()),'simulated_share_north25':float(sim[yy.ravel()>=25].sum()),

            'observed_share_west125':float(obs[xx.ravel()<125].sum()),'simulated_share_west125':float(sim[xx.ravel()<125].sum())})

    pd.DataFrame(baseline).to_csv(_LAYOUT_ROOT / 'results/models/bam_existing_baseline_skill.csv',index=False)

    residual=a['d_full']-a['d_steer']-a['d_loc']-a['d_life']

    g=domain(primary).sort_values(['sid','iso_time']).groupby('sid')

    gap=(g.iso_time.diff()/pd.Timedelta(hours=1)).dropna()

    write_json('results/models/bam_existing_diagnostics.json',{'observed_reconstruction_maxabs':error,

        'observed_tv':tv(d),'single_factor_sum_residual_tv':tv(residual),

        'scenario_keys':a.files,'non_six_hour_segments':int(gap.ne(6).sum()),'n_segments':int(len(gap)),

        'gap_hours_counts':{str(k):int(v) for k,v in gap.value_counts().sort_index().items()},

        'interpretation':'Separate scenario differences are not additive contributions. Norm shares and Shapley projection fractions use different definitions.',

        'source_code_result_mapping':'Reference f_BASE_P1/P2 fields correspond to the f_BASE_EARLY/LATE simulation scenarios.'})



def absolute_metrics(catalogues,primary):

    annual=[]; genesis=[]; recurrence=[]

    collections={'PRIMARY':primary,**{k:v['ts_points'] for k,v in catalogues.items()}}

    for agency,ts in collections.items():

        pts=domain(ts)

        years=range(1966,2026)

        for year in years:

            p=pts.loc[pts.season.eq(year)]

            for low,high in [(0,20),(20,40)]:

                b=p.loc[p.lat.ge(low)&(p.lat.lt(high) if high<40 else p.lat.le(high))]

                annual.append({'agency':agency,'year':year,'latitude_band':f'{low}-{high}',

                    'ts_track_points':len(b),'unique_storms':b.sid.nunique(),

                    'relative_point_share':len(b)/len(p) if len(p) else np.nan})

        # Same inherited sustained-recurvature criterion; native TS eligibility already applied.

        # Agency wind durations are not converted. The threshold below is a dummy pass-through.

        rec_rows=[]

        for sid,g in ts.groupby('sid',sort=False):

            g=g.sort_values('iso_time').copy()

            g['wind']=34.0

            if 'nature' in g:g=g.drop(columns='nature')

            rec,when,lat,lon=recurvature(g,34)

            rec_rows.append({'agency':agency,'sid':sid,'year':int(g.season.iloc[0]),'recurving':int(rec),

                             'recurve_time':when,'recurve_lat':lat,'recurve_lon':lon})

        recurrence.extend(rec_rows)

        first=pts.sort_values(['sid','iso_time']).groupby('sid',sort=False).first().reset_index()

        for period,start,end in [('early',1966,1995),('late',1996,2025)]:

            p=first.loc[first.season.between(start,end)]

            genesis.append({'agency':agency,'period':period,'start':start,'end':end,'storms':len(p),

                'mean_genesis_lat':p.lat.mean(),'mean_genesis_lon':p.lon.mean(),

                'share_south15':p.lat.lt(15).mean(),

                'definition':'first native TS-qualified point inside 100-180E,0-40N; PRIMARY uses LAT/LON+USA wind'})

    ann=pd.DataFrame(annual)

    ann.to_csv(_LAYOUT_ROOT / 'results/tracks/absolute_track_annual.csv',index=False)

    summaries=[]

    for (agency,band),g in ann.groupby(['agency','latitude_band'],sort=False):

        for label,start,end in [('early',1966,1995),('late',1996,2025)]:

            p=g.loc[g.year.between(start,end)]

            summaries.append({'agency':agency,'latitude_band':band,'period':label,'n_years':len(p),

                **{c:float(p[c].mean()) for c in ['ts_track_points','unique_storms','relative_point_share']}})

    pd.DataFrame(summaries).to_csv(_LAYOUT_ROOT / 'results/tracks/absolute_track_periods.csv',index=False)

    pd.DataFrame(genesis).to_csv(_LAYOUT_ROOT / 'results/genesis/genesis_summary.csv',index=False)

    rec=pd.DataFrame(recurrence)

    rec.to_csv(_LAYOUT_ROOT / 'results/tracks/recurvature_storms.csv',index=False)

    result=[]

    for agency,g in rec.groupby('agency',sort=False):

        for label,start,end in [('early',1966,1995),('late',1996,2025)]:

            p=g.loc[g.year.between(start,end)]; r=p.loc[p.recurving.eq(1)]

            result.append({'agency':agency,'period':label,'start':start,'end':end,'eligible_storms':len(p),

                'recurving_storms':len(r),'recurving_per_year':len(r)/(end-start+1),

                'recurving_fraction':len(r)/len(p),'mean_recurve_lat':r.recurve_lat.mean(),'mean_recurve_lon':r.recurve_lon.mean()})

    pd.DataFrame(result).to_csv(_LAYOUT_ROOT / 'results/tracks/recurvature_summary.csv',index=False)

    write_json('results/diagnostics/existing_evidence.json',{'status':'PASS','bam_observed_reproduced':True,

        'agencies':list(collections),'hko_2025_completed_from_official_source':True,

        'recurvature_definition':'inherited analysis.tracks.track_metrics.recurvature; native TS eligibility, no cross-agency wind conversion',

        'limits':['Descriptive metrics, not causal contributions.','HKO includes official 2025 completion; four-agency inference is a separate P4 output.',

                  'A storm can cross multiple latitude bands; unique-storm band counts must not be added.','No claim of an absolute increase without evaluating the relevant count.']})



def main():

    primary=primary_ts()

    bam_report(primary)

    print('BAM observed contrast and scenario diagnostics saved.',flush=True)

    cats=all_catalogues()

    absolute_metrics(cats,primary)

    print('P1 absolute, genesis and recurrence metrics saved for PRIMARY/USA/JMA/CMA/HKO.',flush=True)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
