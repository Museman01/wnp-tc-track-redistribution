'Track baseline diagnostics.'
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import json

import numpy as np

import pandas as pd

import netCDF4 as nc

from core.paths import PROJECT,output_path,write_json,import_source

from analysis.tracks.track_metrics import recurvature



R=6371000.0;DT=21600.0



def climatologies(path):

    with nc.Dataset(path) as ds:

        tv=ds['time'];dates=nc.num2date(tv[:],tv.units,getattr(tv,'calendar','standard'))

        year=np.array([d.year for d in dates]);month=np.array([d.month for d in dates])

        lat=np.asarray(ds['latitude'][:],float);lon=np.asarray(ds['longitude'][:],float)

        out={}

        for label,(a,b) in {'early':(1966,1995),'late':(1996,2025)}.items():

            u=np.empty((12,len(lat),len(lon)));v=np.empty_like(u)

            for m in range(1,13):

                idx=np.flatnonzero((year>=a)&(year<=b)&(month==m))

                u[m-1]=np.asarray(ds['u_steer'][idx],dtype=np.float32).mean(0,dtype=np.float64)

                v[m-1]=np.asarray(ds['v_steer'][idx],dtype=np.float32).mean(0,dtype=np.float64)

            out[label]=(u,v)

    return out,lat,lon



def interp(f,lat,lon,plat,plon):

    """Vectorised form of the BAM bilinear sampler (NaN outside the field)."""

    dy=lat[1]-lat[0];dx=lon[1]-lon[0]

    fy=(plat-lat[0])/dy;fx=(plon-lon[0])/dx

    j=np.floor(fy).astype(int);i=np.floor(fx).astype(int)

    ok=(j>=0)&(i>=0)&(j<len(lat)-1)&(i<len(lon)-1)

    out=np.full(plat.shape,np.nan)

    jj,ii=j[ok],i[ok];y=fy[ok]-jj;x=fx[ok]-ii

    out[ok]=f[jj,ii]*(1-y)*(1-x)+f[jj,ii+1]*(1-y)*x+f[jj+1,ii]*y*(1-x)+f[jj+1,ii+1]*y*x

    return out



def simulate(seeds,u,v,lat,lon,par):

    """Same step order as the archived engine: record, sample genesis-month field, advance."""

    a,ub,vb=par

    plat=seeds.lat.to_numpy(float).copy();plon=seeds.lon.to_numpy(float).copy()

    steps=seeds.steps.to_numpy(int);month=seeds.month.to_numpy(int)-1

    alive=np.ones(len(seeds),bool);records=[]

    for step in range(int(steps.max())):

        alive&=step<steps

        alive&=(plon>=100)&(plon<=180)&(plat>=0)&(plat<=40)

        if not alive.any():break

        k=np.flatnonzero(alive)

        records.append(pd.DataFrame({'k':k,'step':step,'lat':plat[k],'lon':plon[k]}))

        su=np.empty(len(k));sv=np.empty(len(k))

        for m in np.unique(month[k]):

            s=month[k]==m

            su[s]=interp(u[m],lat,lon,plat[k][s],plon[k][s]);sv[s]=interp(v[m],lat,lon,plat[k][s],plon[k][s])

        good=np.isfinite(su)&np.isfinite(sv)

        alive[k[~good]]=False;k=k[good];su=su[good];sv=sv[good]

        cu=a*su+ub;cv=a*sv+vb

        plon[k]+=cu*DT/(R*np.cos(np.deg2rad(plat[k])))*180/np.pi

        plat[k]+=cv*DT/R*180/np.pi

    t=pd.concat(records,ignore_index=True)

    t['sid']=seeds.sid.to_numpy()[t.k];t['season']=seeds.season.to_numpy()[t.k]

    t['iso_time']=seeds.iso_time.to_numpy()[t.k]+pd.to_timedelta(6*t.step,unit='h')

    return t.sort_values(['sid','step'])



def period_field(t,years):

    xe=np.arange(100,182.5,2.5);ye=np.arange(0,42.5,2.5);rows=[]

    for y in years:

        a=t.loc[t.season.eq(y)]

        h=np.histogram2d(a.lon,a.lat,bins=[xe,ye])[0].T

        rows.append(h/h.sum())

    return np.mean(rows,0)



def storm_metrics(t):

    """Storm-level descriptors comparable between observed and synthetic in-domain tracks."""

    rec=[];net=[];useg=[]

    for sid,g in t.groupby('sid',sort=False):

        g=g.sort_values('iso_time')

        x=g[['iso_time','lat','lon']].copy();x['wind']=34.0

        r=recurvature(x,34)[0];rec.append(r)

        net.append(g.lon.iloc[-1]-g.lon.iloc[0])

        if len(g)>1:

            dt=np.diff(g.iso_time.to_numpy())/np.timedelta64(1,'s')

            la=g.lat.to_numpy();lo=g.lon.to_numpy()

            du=np.diff(lo)*np.pi/180*R*np.cos(np.deg2rad((la[1:]+la[:-1])/2))/dt

            ok=(dt>0)&(dt<=43200);useg.extend(du[ok])

    return {'storms':len(rec),'recurving_fraction':float(np.mean(rec)),

            'net_westward_fraction':float(np.mean(np.array(net)<0)),'mean_zonal_motion_ms':float(np.mean(useg))}



def bam_diagnostics():

    bam=import_source('baseline_motion_reader','analysis/models/counterfactual_steering.py')

    cfg=bam.load_config();tracks,_=bam.load_tracks(cfg)

    g=tracks.groupby('sid')

    seeds=g.first().reset_index()[['sid','season','lat','lon','iso_time']]

    seeds['month']=seeds.iso_time.dt.month;seeds['steps']=g.size().reindex(seeds.sid).to_numpy()

    clim,lat,lon=climatologies(PROJECT/'data/interim/steering.nc')

    archived=np.load(_LAYOUT_ROOT / 'results/models/counterfactual/fields.npz')

    old=json.loads((_LAYOUT_ROOT / 'results/models/counterfactual/criteria.json').read_text())['calibration']

    sens=pd.read_csv(_LAYOUT_ROOT / 'results/models/bam_sensitivity_scenarios.csv')

    act=sens.loc[sens.case.eq('850-300_actual_dt_max12h')].iloc[0]

    rounds=np.load(_LAYOUT_ROOT / 'results/models/bam_sensitivity_annual.npz')

    cases={'original_calibration':((old['alpha'],old['u_beta'],old['v_beta']),{'early':archived['f_BASE_P1'],'late':archived['f_BASE_P2']}),

           '850-300_actual_dt_max12h':((act.alpha,act.u_intercept_ms,act.v_intercept_ms),

                {'early':rounds['850-300_actual_dt_max12h_early'].mean(0).reshape(16,32),'late':rounds['850-300_actual_dt_max12h_late'].mean(0).reshape(16,32)})}

    rows=[];checks={}

    periods={'early':(1966,1995),'late':(1996,2025)}

    for label,(a,b) in periods.items():

        obs=tracks.loc[tracks.season.between(a,b)]

        f=period_field(obs,range(a,b+1))

        rows.append({'case':'observed','period':label,'share_north20':float(f[8:].sum()),**storm_metrics(obs)})

    for case,(par,ref) in cases.items():

        for label,(a,b) in periods.items():

            s=seeds.loc[seeds.season.between(a,b)]

            t=simulate(s,*clim[label],lat,lon,par)

            f=period_field(t,range(a,b+1))

            err=float(np.abs(f-ref[label]).max());checks[f'{case}_{label}']=err

            if err>1e-10:raise AssertionError(f'Regenerated {case} {label} field differs: {err}')

            rows.append({'case':case,'period':label,'share_north20':float(f[8:].sum()),**storm_metrics(t)})

        print('BAM storm diagnostics',case,flush=True)

    pd.DataFrame(rows).to_csv(_LAYOUT_ROOT / 'results/models/bam_baseline_track_diagnostics.csv',index=False)

    # Latitude-band share for every round-1 variant, using the Table S8 20N division.

    band=[]

    obs=rounds['observed']

    for period,sl in [('early',slice(0,30)),('late',slice(30,60))]:

        band.append({'case':'observed','period':period,'share_north20':float(obs[sl].mean(0).reshape(16,32)[8:].sum())})

    for key in rounds.files:

        for period in ['early','late']:

            if key.endswith('_'+period):

                band.append({'case':key[:-len(period)-1],'period':period,'share_north20':float(rounds[key].mean(0).reshape(16,32)[8:].sum())})

    pd.DataFrame(band).to_csv(_LAYOUT_ROOT / 'results/models/bam_baseline_band_shares.csv',index=False)

    return checks



def landfall_bands():

    e=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/hko_landfall_events.csv',parse_dates=['time'])

    e=e.loc[e.pre_crossing_native_ts&e.season.between(1966,2025)]

    rows=[]

    for kind in ['all_events','first_landfall']:

        x=e if kind=='all_events' else e.sort_values('time').drop_duplicates(['agency','sid'])

        for agency,g in x.groupby('agency'):

            for band,mask in [('0-20N',g.lat.lt(20)),('20N_north',g.lat.ge(20))]:

                b=g.loc[mask];annual=b.groupby('season').size().reindex(range(1966,2026),fill_value=0)

                rows.append({'agency':agency,'definition':kind,'latitude_band':band,

                    'early_per_year':float(annual.loc[1966:1995].mean()),'late_per_year':float(annual.loc[1996:2025].mean()),

                    'difference_per_year':float(annual.loc[1996:2025].mean()-annual.loc[1966:1995].mean())})

    pd.DataFrame(rows).to_csv(_LAYOUT_ROOT / 'results/landfall/landfall_latitude_band_counts.csv',index=False)



def ridge_context():

    w=pd.read_csv(PROJECT/'data/processed/subtropical_high_indices.csv')

    rows=[]

    for c in ['ridge_line','west_ridge_point','wpsh_area','wpsh_intensity']:

        e=w.loc[w.season.between(1966,1995),c];l=w.loc[w.season.between(1996,2025),c]

        rows.append({'metric':c,'early_mean':float(e.mean()),'early_defined_years':int(e.notna().sum()),

            'late_mean':float(l.mean()),'late_defined_years':int(l.notna().sum()),'difference':float(l.mean()-e.mean())})

    pd.DataFrame(rows).to_csv(_LAYOUT_ROOT / 'results/circulation/subtropical_high_period_means.csv',index=False)



def main():

    checks=bam_diagnostics();landfall_bands();ridge_context()

    write_json('results/diagnostics/track_diagnostics.json',{'status':'PASS',

        'bam_field_regeneration_maxabs':checks,

        'bam_storm_metrics':'In-domain TS-qualified observed tracks versus regenerated baseline trajectories; sustained recurvature uses analysis.tracks.track_metrics.recurvature on in-domain positions; zonal motion uses 0<gap<=12 h segments.',

        'landfall_bands':'Four-agency exact-crossing events with a 20N latitude division; descriptive counts only.',

        'ridge':'data/processed/subtropical_high_indices.csv; ridge_line is the wind-defined U500 zero-crossing index; fixed-contour metrics undefined when the 5880-gpm contour is absent.'})

    print('Round-2 diagnostics saved',flush=True)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
