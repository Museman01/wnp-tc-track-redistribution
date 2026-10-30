'Hko figure series.'
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import sys,json

import numpy as np,pandas as pd

from shapely.geometry import LineString,MultiLineString,Point

from shapely.ops import transform

from core.paths import PROJECT,output_path,write_json



import analysis.landfall.native_landfall_validation as LL                 # noqa: E402

import analysis.tracks.common_storm_validation as CS               # noqa: E402

import analysis.tracks.cutpoint_validation as CP              # noqa: E402

import analysis.landfall.coast_grouping_validation as CG        # noqa: E402

from analysis.shared.block_resampling import block_bootstrap_many, block_permutation_many  # noqa: E402

import analysis.shared.record_composition as RC

from analysis.shared.agency_catalogues import all_catalogues                            # noqa: E402



CFG={'n_permutations':9999,'random_seed':202406,'n_bootstrap':4999}

Z=_LAYOUT_ROOT / 'data/reference'

OUT=_LAYOUT_ROOT / 'results/agency_comparison/figure_series';YEARS=np.arange(1966,2026)

checks={};info_extra={}



def close(a,b,tol=1e-9):return bool(np.allclose(np.asarray(a,float),np.asarray(b,float),atol=tol,equal_nan=True))



def hko_events():

    e=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/hko_landfall_events.csv')

    return e.loc[e.pre_crossing_native_ts.astype(bool)].copy()



def formal_events():

    e=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/landfall_events_exact.csv')

    return e.loc[e.pre_crossing_native_ts.astype(bool)].copy()



def frequency_and_lmi(cats):

    core=pd.read_csv(Z/'agency_comparison/03_common_storms/core_crossagency_annual.csv')

    # Reproduce the three-agency series with the formal catalog builder.

    common=set.intersection(*(set(cats[a]['lmi']['sid']) for a in ['USA','TOKYO','CMA']))

    for a in ['USA','TOKYO','CMA']:

        f=cats[a]['frequency'].set_index('season')['n_tc'].reindex(YEARS,fill_value=0)

        ref=core[core.agency.eq(a)].set_index('season')

        checks[f'frequency_{a}']=close(f.values,ref.loc[YEARS,'n_tc'].values)

        lmi=cats[a]['lmi'];lmi=lmi[lmi.season.between(1982,2025)&lmi.sid.isin(common)].groupby('season')['lmi_lat'].mean().reindex(range(1982,2026))

        checks[f'lmi_common_{a}']=close(lmi.values,ref.loc[1982:2025,'mean_lmi_lat_common'].values)

    h=cats['HKO']['ts_points']

    freq=h.groupby('season')['sid'].nunique().reindex(YEARS,fill_value=0)

    w=h.loc[h.wind.notna()];lmi=w.loc[w.groupby('sid')['wind'].idxmax()]

    lmi=lmi[lmi.season.between(1982,2025)&lmi.sid.isin(common)]

    m=lmi.groupby('season')['lat'].mean().reindex(range(1982,2026))

    pd.DataFrame({'season':YEARS,'n_tc':freq.values}).to_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/frequency.csv',index=False)

    pd.DataFrame({'season':range(1982,2026),'mean_lmi_lat_common':m.values,'n_storms':lmi.groupby('season').size().reindex(range(1982,2026),fill_value=0).values}).to_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/lmi_common.csv',index=False)

    return dict(hko_lmi_common_storms=int(len(lmi)),three_agency_common_1982_2025=int(sum(1 for s in common if 1982<=int(s[:4])<=2025)))



def index(cats):

    ref=pd.read_csv(Z/'environment_validation/supplemental_validation/results/wnp_tc_redistribution_index_annual.csv')

    usa=RC.load_agency_tracks()['USA']

    point,*_=RC.annual_fields(usa,YEARS,2.5);_,oos,_=RC.projection_scores(point,YEARS)

    r=ref[ref.agency.eq('USA')&ref.weighting.eq('track_point')].sort_values('year')

    checks['index_USA']=close(oos,r.index_oos.values,1e-9)

    h=cats['HKO']['ts_points'];h=h[h.season.between(1966,2025)&h.lon.between(100,180)&h.lat.between(0,40)]

    point,*_=RC.annual_fields(h,YEARS,2.5);_,oos,_=RC.projection_scores(point,YEARS)

    pd.DataFrame({'year':YEARS,'agency':'HKO','weighting':'track_point','index_oos':oos}).to_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/index.csv',index=False)

    return dict(hko_index_late_minus_early=float(oos[30:].mean()-oos[:30].mean()))



def landfall():

    ev=pd.concat([formal_events(),hko_events().loc[lambda d:d.agency.eq('HKO')]],ignore_index=True)

    rows=[]

    for ag,g in ev.groupby('agency'):

        for dfn,fr in LL.definitions(g.assign(coast=g.get('coast','Other'))).items():

            if dfn=='unique_storm_coast':continue

            for r in LL.annual_metrics(fr,YEARS):rows.append({'agency':ag,'definition':dfn,**r})

    ann=pd.DataFrame(rows)

    out=LL.summarize(ann,CFG);out=out[(out.start==1966)&(out.end==2025)&(out.metric=='mean_lat')]

    ref=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/landfall_latitude_summary.csv');ref=ref[(ref.start==1966)&(ref.end==2025)&(ref.metric=='mean_lat')]

    m=out.merge(ref,on=['agency','definition'],suffixes=('','_ref'))

    checks['landfall_three_agency']=close(m[['period_difference','period_ci_low','period_ci_high','period_block_p']].values,

                                          m[['period_difference_ref','period_ci_low_ref','period_ci_high_ref','period_block_p_ref']].values)

    out[['agency','definition','period_difference','period_ci_low','period_ci_high','period_block_p']].to_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/landfall.csv',index=False)

    return ev



def common_storms(ev,cats,prefilter):

    common=set.intersection(*(set(cats[a]['ts_points'].loc[cats[a]['ts_points'].season.between(1966,2025),'sid']) for a in ['USA','TOKYO','CMA']))

    if not prefilter:ev=pd.concat([pd.read_csv(_LAYOUT_ROOT / 'results/landfall/landfall_events_exact.csv'),hko_events().loc[lambda d:d.agency.eq('HKO')]],ignore_index=True)

    ref=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/common_storm_landfall_summary.csv');ref=ref[ref.metric.eq('mean_landfall_latitude')]

    e=ev.copy();e['coast']=e.get('coast','Other');e['coast']=e['coast'].fillna('Other')

    f=e.loc[e.sid.isin(common)].copy();f['time']=pd.to_datetime(f.time)

    first=f.sort_values('time').drop_duplicates(['agency','sid'],keep='first')

    early=np.flatnonzero(YEARS<=1995);late=np.flatnonzero(YEARS>=1996);rows=[]

    for ag in ['USA','TOKYO','CMA','HKO']:

        for dfn,fr in [('all_events',f),('first_landfall',first)]:

            part=fr.loc[fr.agency.eq(ag)]

            lat=part.groupby('season')['lat'].mean().reindex(YEARS)

            named=part.loc[part.coast.isin(CS.NORTH|CS.SOUTH)]

            n=named.coast.isin(CS.NORTH).groupby(named.season).sum().reindex(YEARS,fill_value=0)

            s=(~named.coast.isin(CS.NORTH)).groupby(named.season).sum().reindex(YEARS,fill_value=0)

            share=(n/(n+s)).where((n+s)>0)

            matrix=np.column_stack([lat.values,share.values]).astype(float)

            ch,p=block_permutation_many(matrix,early,late,block=3,nperm=CFG['n_permutations'],seed=CFG['random_seed'])

            lo,hi=block_bootstrap_many(matrix[early],matrix[late],block=3,nboot=CFG['n_bootstrap'],seed=CFG['random_seed'])

            rows.append(dict(agency=ag,definition=dfn,late_minus_early=ch[0],ci_low=lo[0],ci_high=hi[0],block_permutation_p=p[0]))

    out=pd.DataFrame(rows)

    m=out.merge(ref,on=['agency','definition'],suffixes=('','_ref'))

    checks[f'common_three_agency_prefilter_{prefilter}']=close(m[['late_minus_early','ci_low','ci_high','block_permutation_p']].values,

                                        m[['late_minus_early_ref','ci_low_ref','ci_high_ref','block_permutation_p_ref']].values)

    if checks[f'common_three_agency_prefilter_{prefilter}']:out.to_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/common.csv',index=False)

    return checks[f'common_three_agency_prefilter_{prefilter}']



def classify(ev):

    """Nearest named segment and distance (km), local equirectangular metric, vectorized."""

    seg=json.loads(_LAYOUT_ROOT / 'results/tracks/fig01_coast_segments.json'.read_text())['segments']

    pieces={k:np.array([[a[0],a[1],b[0],b[1]] for l in v for a,b in zip(l[:-1],l[1:])]) for k,v in seg.items()}

    lon=ev.lon.to_numpy(float);lat=ev.lat.to_numpy(float);keys=list(pieces);best=np.full((len(lon),len(keys)),np.inf)

    for j,k in enumerate(keys):

        P=pieces[k]

        for i0 in range(0,len(lon),100):

            lo=lon[i0:i0+100,None];la=lat[i0:i0+100,None];c=np.cos(np.radians(la))

            x1=(P[None,:,0]-lo)*111.32*c;y1=(P[None,:,1]-la)*110.57;x2=(P[None,:,2]-lo)*111.32*c;y2=(P[None,:,3]-la)*110.57

            dx=x2-x1;dy=y2-y1;L=dx*dx+dy*dy

            t=np.clip(np.where(L>0,-(x1*dx+y1*dy)/np.where(L>0,L,1),0),0,1)

            d=np.hypot(x1+t*dx,y1+t*dy).min(axis=1);best[i0:i0+100,j]=d

    idx=best.argmin(axis=1)

    return np.array(keys)[idx],best[np.arange(len(lon)),idx]



def coast(ev):

    from analysis.shared.block_resampling import block_permutation_difference

    f=ev.loc[ev.agency.isin(['USA','TOKYO','CMA'])].copy()

    near,dist=classify(f)

    mine=np.where(dist<=50,near,'Other');auth=np.where(f.nearest_distance_km.values<=50,f.nearest_segment.values,'Other')

    agree=float((mine==auth).mean());checks['coast_assignment_agreement_formal_events']=agree>=0.999

    info_extra['coast_assignment_agreement']=f'{int((mine==auth).sum())}/{len(mine)}'

    h=ev.loc[ev.agency.eq('HKO')].copy();near,dist=classify(h);h['nearest_segment']=near;h['nearest_distance_km']=dist

    cols=['agency','sid','season','time','lat','lon','nearest_segment','nearest_distance_km']

    allev=pd.concat([f[cols],h[cols]],ignore_index=True);allev['time']=pd.to_datetime(allev.time);allev['wind']=np.nan

    allev['coast']=np.where(allev.nearest_distance_km.astype(float)<=50,allev.nearest_segment,'Other')

    early,late=np.arange(30),np.arange(30,60);rows=[]

    for rule in ['all_events','first_any','first_named']:

        frame=CG.select_events(allev,rule)

        for ag in ['USA','TOKYO','CMA','HKO']:

            share=CG.annual_share(frame,ag,CG.BASE_NORTH,CG.BASE_SOUTH,'named_only')

            eff,p=block_permutation_difference(share,early,late,block=3,nperm=CFG['n_permutations'],seed=CFG['random_seed'])

            rows.append(dict(agency=ag,event_rule=rule,change_percentage_points=eff*100,block_permutation_p=p))

    res=pd.DataFrame(rows)

    ref=pd.read_csv(Z/'agency_comparison/06_landfall_grouping/coast_grouping_sensitivity.csv')

    ref=ref[ref.threshold_km.eq(50.0)&ref.taiwan_rule.eq('taiwan_north')&ref.denominator.eq('named_only')]

    m=res.merge(ref,on=['agency','event_rule'],suffixes=('','_ref'))

    checks['coast_three_agency']=close(m[['change_percentage_points','block_permutation_p']].values,m[['change_percentage_points_ref','block_permutation_p_ref']].values,1e-9)

    res.to_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/coast.csv',index=False)

    h[cols].to_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/hko_event_coast_classification.csv',index=False)

    return h



def cutpoint(ev,hko_classified):

    ref=pd.read_csv(Z/'agency_comparison/02_cutpoint_sensitivity/cutpoint_sensitivity.csv')

    ref=ref[ref.scheme.eq('full_record')&ref.metric.eq('first_landfall_latitude_difference')]

    e=ev.copy();e['coast']=e.get('coast','Other');e['coast']=e['coast'].fillna('Other')

    land=CP.landfall_annual(e,YEARS);rows=[]

    for cut in range(1991,2002):

        early=np.arange(1966,cut);late=np.arange(cut,2026)

        for ag in ['USA','TOKYO','CMA','HKO']:

            lf=land.loc[land.agency.eq(ag)]

            eff,p=CP.scalar_test(lf.season,lf.first_landfall_mean_latitude,early,late,CFG)

            rows.append(dict(cutpoint=cut,agency=ag,effect=eff,p_value=p))

    out=pd.DataFrame(rows);m=out.merge(ref,on=['cutpoint','agency'],suffixes=('','_ref'))

    checks['cutpoint_three_agency']=close(m[['effect','p_value']].values,m[['effect_ref','p_value_ref']].values)

    out.to_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/cutpoint.csv',index=False)



def main(stage):

    OUT.mkdir(parents=True,exist_ok=True);qa=output_path('results/diagnostics/hko_figure_series.json')

    state=json.loads(qa.read_text()) if qa.exists() else {}

    checks.update(state.get('reproduction_checks',{}));info=state.get('info',{})

    cats=all_catalogues() if stage in ('basic','common','decomp') else None

    if stage=='basic':info.update(frequency_and_lmi(cats));info.update(index(cats))

    ev=landfall() if stage in ('landfall','common','coast','cutpoint') else None

    if stage=='common':

        if not common_storms(ev,cats,True):common_storms(ev,cats,False)

    if stage=='coast':coast(ev);info.update(info_extra)

    if stage=='cutpoint':

        if not cutpoint(ev,True):cutpoint(ev,False)

    variant=lambda k:k.startswith('common_') or k.startswith('cutpoint_')

    done={'basic':'index_USA' in checks,'landfall':'landfall_three_agency' in checks,'common':any(k.startswith('common_') for k in checks),

          'coast':'coast_three_agency' in checks,'cutpoint':any(k.startswith('cutpoint_') for k in checks)}

    ok=all(done.values()) and all(v for k,v in checks.items() if not variant(k)) and any(v for k,v in checks.items() if k.startswith('common_')) and any(v for k,v in checks.items() if k.startswith('cutpoint_'))

    write_json('results/diagnostics/hko_figure_series.json',dict(status='PASS' if ok else ('INCOMPLETE' if not all(done.values()) else 'CHECK_FAILED'),

        reproduction_checks=checks,settings=CFG,info=info,

        note='The same statistical definitions are applied to the three reference agency records and HKO.'))

    print(stage,{k:v for k,v in checks.items()},info)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main(sys.argv[1])
