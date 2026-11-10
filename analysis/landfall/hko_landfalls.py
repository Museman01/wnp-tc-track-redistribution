"""Four-agency exact sea-to-land geometry on native tracks, with explicit method reconciliation."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import numpy as np

import pandas as pd

import geopandas as gpd

import shapely

from shapely.geometry import LineString,Point

from core.paths import PROJECT,output_path

from analysis.shared.agency_catalogues import all_catalogues

from analysis.shared.inference_stats import bh_fdr

from analysis.tracks.multiagency_sensitivity import _detect_landfalls

from core.landfall_gshhg import build_land_mask

import yaml

import argparse



def points(g):

    if g.is_empty:return []

    if g.geom_type=='Point':return [g]

    if g.geom_type in ['LineString','LinearRing']:return [Point(g.coords[0]),Point(g.coords[-1])]

    return [p for child in getattr(g,'geoms',[]) for p in points(child)]



def crossings(line,geoms):

    intersections=shapely.intersection(shapely.boundary(geoms),line)

    fractions=sorted(line.project(p,normalized=True) for g in intersections for p in points(g))

    unique=[]

    for f in fractions:

        if not unique or abs(f-unique[-1])>1e-8:unique.append(float(np.clip(f,0,1)))

    if not unique:return []

    bounds=[0.,*unique,1.]

    state=[bool(np.any(shapely.covers(geoms,line.interpolate((a+b)/2 if b-a>1e-10 else min(1,a+1e-8),normalized=True)))) for a,b in zip(bounds[:-1],bounds[1:])]

    return [f for i,f in enumerate(unique) if not state[i] and state[i+1]]



def detect(cat,geoms,tree,mask_info,cfg,start=1966,end=2025):

    t=cat['tracks'].loc[cat['tracks'].season.between(start,end)].sort_values(['sid','iso_time']).copy()

    ts=set(zip(cat['ts_points'].sid,cat['ts_points'].iso_time))

    for c in ['lat','lon','iso_time','sid']:t['next_'+c]=t.groupby('sid')[c].shift(-1)

    t['gap']=(t.next_iso_time-t.iso_time)/pd.Timedelta(hours=1)

    # Domain membership is evaluated at the crossing by the raster screen.

    # A valid coastal crossing can have one endpoint outside the domain.

    t=t.loc[t.sid.eq(t.next_sid)&t.gap.gt(0)&t.gap.le(12)&t[['lat','lon','next_lat','next_lon']].notna().all(axis=1)].reset_index(drop=True)

    candidates=_detect_landfalls(cat['tracks'].loc[cat['tracks'].season.between(start,end)],cat['agency'],cfg,mask_info)

    print('Raster candidates',cat['agency'],len(candidates),flush=True)

    events=[];groups={sid:g for sid,g in t.groupby('sid',sort=False)}

    for candidate in candidates.itertuples():

        group=groups.get(candidate.sid)

        if group is None:continue

        bracket=group.loc[group.iso_time.le(candidate.time)&group.next_iso_time.ge(candidate.time)]

        if bracket.empty:continue

        r=bracket.iloc[0];line=LineString([(r.lon,r.lat),(r.next_lon,r.next_lat)])

        ids=tree.query(line)

        fractions=crossings(line,geoms[np.asarray(ids)]) if len(ids) else []

        if not fractions:continue

        approximate=(candidate.time-r.iso_time)/(r.next_iso_time-r.iso_time)

        f=min(fractions,key=lambda value:abs(value-approximate))

        when=r.iso_time+(r.next_iso_time-r.iso_time)*f

        events.append({'agency':cat['agency'],'sid':r.sid,'season':int(r.season),'time':when,

            'lat':float(r.lat+f*(r.next_lat-r.lat)),'lon':float(r.lon+f*(r.next_lon-r.lon)),

            'pre_crossing_native_ts':(r.sid,r.iso_time) in ts})

    e=pd.DataFrame(events).sort_values(['sid','time']).drop_duplicates(['sid','time'])

    keep=[]

    for _,g in e.groupby('sid',sort=False):

        last=None

        for idx,r in g.iterrows():

            if last is None or (r.time-last)/pd.Timedelta(hours=1)>=12:keep.append(idx);last=r.time

    return e.loc[keep].loc[lambda x:x.pre_crossing_native_ts].copy()



def inference(annual):

    vals=np.asarray(annual,float);obs=vals[30:].mean()-vals[:30].mean()

    n=len(vals);nlate=n-30

    rng=np.random.default_rng(202406);blocks=[np.arange(i,min(i+3,n)) for i in range(0,n,3)]

    exceed=0

    for _ in range(9999):

        order=np.concatenate([blocks[i] for i in rng.permutation(len(blocks))])

        exceed+=abs(vals[order[30:]].mean()-vals[order[:30]].mean())>=abs(obs)

    boot=[]

    for _ in range(1999):

        ix=np.concatenate([(start+np.arange(3))%30 for start in rng.integers(0,30,10)])

        jx=np.concatenate([(start+np.arange(3))%nlate for start in rng.integers(0,nlate,int(np.ceil(nlate/3)))])[:nlate]

        boot.append(vals[30:][jx].mean()-vals[:30][ix].mean())

    low,high=np.quantile(boot,[.025,.975])

    return obs,(exceed+1)/10000,low,high



def main():

    parser=argparse.ArgumentParser();parser.add_argument('--recompute-usa',action='store_true');parser.add_argument('--extend-2025',action='store_true');parser.add_argument('--output-dir',default='results/landfall');args=parser.parse_args()
    out=(PROJECT / args.output_dir).resolve()
    if not out.is_relative_to(PROJECT.resolve()):raise ValueError('Output directory must remain within the project')
    names=['hko_landfall_events.csv','hko_landfall_summary.csv','hko_landfall_annual.csv','hko_landfall_method_reconciliation.csv']
    cache=_LAYOUT_ROOT / 'results/landfall/hko_landfall_events.csv'
    if (args.recompute_usa or args.extend_2025) and not cache.is_file():raise FileNotFoundError(cache)
    required=names if args.recompute_usa or args.extend_2025 or not cache.exists() else names[1:]
    if not (PROJECT/'data/derived/land_mask.npz').exists():required=required+['land_mask.npz']
    if any((out/name).exists() for name in required):raise FileExistsError('Refusing to overwrite existing landfall results')
    out.mkdir(parents=True,exist_ok=True)

    cats=all_catalogues()

    if cache.exists() and not args.recompute_usa and not args.extend_2025:events=pd.read_csv(cache,parse_dates=['time'])

    else:

        coast=gpd.read_file(PROJECT/'data/raw/GSHHG/GSHHS_shp/h/GSHHS_h_L1.shp',bbox=(80,0,180,65))

        geoms=coast.geometry.to_numpy();tree=shapely.STRtree(geoms)

        cfg=yaml.safe_load((PROJECT/'config/config.yaml').read_text(encoding='utf-8'))

        cfg['gshhg_path']=str(PROJECT/cfg['gshhg_path'])

        maskpath=output_path('data/derived/land_mask.npz')

        if maskpath.exists():

            m=np.load(maskpath);mask_info=(m['mask'],80.,65.,.02)

        else:

            mask_info=build_land_mask(cfg);np.savez_compressed(out/'land_mask.npz',mask=mask_info[0])

        print('Land mask and vector coast ready',flush=True)

        frames=[]

        if args.extend_2025:

            prior=pd.read_csv(cache,parse_dates=['time']);frames.append(prior.loc[prior.season.le(2024)])

            formal=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/landfall_events_exact.csv',parse_dates=['time'])

            formal=formal.loc[formal.season.eq(2025)&formal.pre_crossing_native_ts]

            frames.append(formal[['agency','sid','season','time','lat','lon','pre_crossing_native_ts']])

        elif args.recompute_usa and cache.exists():

            prior=pd.read_csv(cache,parse_dates=['time'])

            frames.append(prior.loc[~prior.agency.eq('USA')])

        for agency,cat in cats.items():

            if args.extend_2025 and agency!='HKO':continue

            if args.recompute_usa and agency!='USA':continue

            e=detect(cat,geoms,tree,mask_info,cfg,start=2025 if args.extend_2025 else 1966);frames.append(e);print('Exact-vector landfalls',agency,len(e),flush=True)

        events=pd.concat(frames,ignore_index=True);events.to_csv(out/'hko_landfall_events.csv',index=False)

    common=set.intersection(*(set(c['eligible_sids']) for c in cats.values()))

    rows=[];yearrows=[]

    for sample in ['all_native','four_agency_common_storms']:

        for agency in cats:

            e=events.loc[events.agency.eq(agency)]

            if sample!='all_native':e=e.loc[e.sid.isin(common)]

            for kind,p in [('all_events',e),('first_landfall',e.sort_values('time').drop_duplicates('sid'))]:

                annual=p.groupby('season').lat.mean().reindex(range(1966,2026))

                if annual.isna().any():raise ValueError(f'Missing annual landfalls: {sample},{agency},{kind}')

                delta,pvalue,lo,hi=inference(annual.to_numpy())

                rows.append({'sample':sample,'agency':agency,'definition':kind,'n_events':len(p),

                    'equal_year_lat_early':annual.iloc[:30].mean(),'equal_year_lat_late':annual.iloc[30:].mean(),

                    'difference_deg':delta,'block_permutation_p':pvalue,'q025_deg':lo,'q975_deg':hi})

                for year,value in annual.items():yearrows.append({'sample':sample,'agency':agency,'definition':kind,'year':year,'mean_lat':value})

    r=pd.DataFrame(rows)

    for _,ids in r.groupby(['sample','definition']).groups.items():r.loc[ids,'q_bh_four_agency_family']=bh_fdr(r.loc[ids,'block_permutation_p'])

    r.to_csv(out / 'hko_landfall_summary.csv',index=False)

    pd.DataFrame(yearrows).to_csv(out / 'hko_landfall_annual.csv',index=False)

    old=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/landfall_events_exact.csv',parse_dates=['time'])

    old=old.loc[old.season.between(1966,2025)&old.pre_crossing_native_ts]

    reconcile=[]

    for agency in ['USA','TOKYO','CMA']:

        a=old.loc[old.agency.eq(agency)];b=events.loc[events.agency.eq(agency)]

        reconcile.append({'agency':agency,'formal_events_same_window':len(a),'direct_vector_events':len(b),

            'difference':len(b)-len(a)})

    pd.DataFrame(reconcile).to_csv(out / 'hko_landfall_method_reconciliation.csv',index=False)

if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()

