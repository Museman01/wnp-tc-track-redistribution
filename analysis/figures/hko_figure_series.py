'Hko figure series.'
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import sys,json

import numpy as np,pandas as pd

from pathlib import Path



import analysis.landfall.native_landfall_validation as LL                 # noqa: E402

import analysis.tracks.common_storm_validation as CS               # noqa: E402

import analysis.tracks.cutpoint_validation as CP              # noqa: E402

import analysis.landfall.coast_grouping_validation as CG        # noqa: E402

from analysis.shared.block_resampling import block_bootstrap_many, block_permutation_many  # noqa: E402

import analysis.shared.record_composition as RC

from analysis.shared.agency_catalogues import all_catalogues                            # noqa: E402



CFG={'n_permutations':9999,'random_seed':202406,'n_bootstrap':4999}

YEARS=np.arange(1966,2026)



def hko_events():

    e=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/hko_landfall_events.csv')

    return e.loc[e.pre_crossing_native_ts.astype(bool)].copy()



def formal_events():

    e=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/landfall_events_exact.csv')

    return e.loc[e.pre_crossing_native_ts.astype(bool)].copy()



def frequency_and_lmi(cats,destination):

    common=set.intersection(*(set(cats[a]['lmi']['sid']) for a in ['USA','TOKYO','CMA']))

    h=cats['HKO']['ts_points']

    freq=h.groupby('season')['sid'].nunique().reindex(YEARS,fill_value=0)

    w=h.loc[h.wind.notna()];lmi=w.loc[w.groupby('sid')['wind'].idxmax()]

    lmi=lmi[lmi.season.between(1982,2025)&lmi.sid.isin(common)]

    m=lmi.groupby('season')['lat'].mean().reindex(range(1982,2026))

    pd.DataFrame({'season':YEARS,'n_tc':freq.values}).to_csv(destination / 'frequency.csv',index=False)

    pd.DataFrame({'season':range(1982,2026),'mean_lmi_lat_common':m.values,'n_storms':lmi.groupby('season').size().reindex(range(1982,2026),fill_value=0).values}).to_csv(destination / 'lmi_common.csv',index=False)

    return dict(hko_lmi_common_storms=int(len(lmi)),three_agency_common_1982_2025=int(sum(1 for s in common if 1982<=int(s[:4])<=2025)))



def index(cats,destination):

    h=cats['HKO']['ts_points'];h=h[h.season.between(1966,2025)&h.lon.between(100,180)&h.lat.between(0,40)]

    point,*_=RC.annual_fields(h,YEARS,2.5);_,oos,_=RC.projection_scores(point,YEARS)

    pd.DataFrame({'year':YEARS,'agency':'HKO','weighting':'track_point','index_oos':oos}).to_csv(destination / 'index.csv',index=False)

    return dict(hko_index_late_minus_early=float(oos[30:].mean()-oos[:30].mean()))



def landfall(destination):

    ev=pd.concat([formal_events(),hko_events().loc[lambda d:d.agency.eq('HKO')]],ignore_index=True)

    rows=[]

    for ag,g in ev.groupby('agency'):

        for dfn,fr in LL.definitions(g.assign(coast=g.get('coast','Other'))).items():

            if dfn=='unique_storm_coast':continue

            for r in LL.annual_metrics(fr,YEARS):rows.append({'agency':ag,'definition':dfn,**r})

    ann=pd.DataFrame(rows)

    out=LL.summarize(ann,CFG);out=out[(out.start==1966)&(out.end==2025)&(out.metric=='mean_lat')]

    out[['agency','definition','period_difference','period_ci_low','period_ci_high','period_block_p']].to_csv(destination / 'landfall.csv',index=False)

    return ev



def common_storms(ev,cats,destination):

    common=set.intersection(*(set(cats[a]['ts_points'].loc[cats[a]['ts_points'].season.between(1966,2025),'sid']) for a in ['USA','TOKYO','CMA']))

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

    out.to_csv(destination / 'common.csv',index=False)



def classify(ev):

    """Nearest named segment and distance (km), local equirectangular metric, vectorized."""

    seg=json.loads((_LAYOUT_ROOT / 'results/tracks/fig01_coast_segments.json').read_text())['segments']

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



def coast(ev,destination):

    from analysis.shared.block_resampling import block_permutation_difference

    f=ev.loc[ev.agency.isin(['USA','TOKYO','CMA'])].copy()

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

    res.to_csv(destination / 'coast.csv',index=False)

    h[cols].to_csv(destination / 'hko_event_coast_classification.csv',index=False)

    return h



def cutpoint(ev,destination):

    e=ev.copy();e['coast']=e.get('coast','Other');e['coast']=e['coast'].fillna('Other')

    land=CP.landfall_annual(e,YEARS);rows=[]

    for cut in range(1991,2002):

        early=np.arange(1966,cut);late=np.arange(cut,2026)

        for ag in ['USA','TOKYO','CMA','HKO']:

            lf=land.loc[land.agency.eq(ag)]

            eff,p=CP.scalar_test(lf.season,lf.first_landfall_mean_latitude,early,late,CFG)

            rows.append(dict(cutpoint=cut,agency=ag,effect=eff,p_value=p))

    out=pd.DataFrame(rows)

    out.to_csv(destination / 'cutpoint.csv',index=False)



def main(stage, output):
    destination=(_LAYOUT_ROOT / output).resolve()
    if not destination.is_relative_to(_LAYOUT_ROOT.resolve()): raise ValueError('Output must remain within the project')
    if destination.exists(): raise FileExistsError('Provide a new output directory')
    if stage not in ('basic','landfall','common','coast','cutpoint'): raise ValueError('Unknown stage')
    destination.mkdir(parents=True)
    info={}
    cats=all_catalogues() if stage in ('basic','common') else None
    if stage=='basic':
        info.update(frequency_and_lmi(cats,destination));info.update(index(cats,destination))
    ev=landfall(destination) if stage in ('landfall','common','coast','cutpoint') else None
    if stage=='common': common_storms(ev,cats,destination)
    if stage=='coast': coast(ev,destination)
    if stage=='cutpoint': cutpoint(ev,destination)
    print(stage,info)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['basic','landfall','common','coast','cutpoint'])
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--output-dir',required=True)
    args=parser.parse_args()
    if not args.execute: parser.error('add --execute to calculate figure data')
    main(args.stage,args.output_dir)
