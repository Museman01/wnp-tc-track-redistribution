"""Six-hourly, storm-centre wind-driven trajectory baseline diagnostics.

The supplied fitted parameters describe 1966-2025. Integration requires external
six-hourly ERA5 fields and includes storm circulation. It is not a
vortex-removed environmental-steering attribution experiment.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import xarray as xr
from numba import njit
from core.paths import PROJECT
from analysis.models.bam_sensitivity import interp,fit,skill
from analysis.tracks.track_baseline_diagnostics import storm_metrics
from analysis.shared.agency_catalogues import annual_fields as repo_annual_fields
YEARS=np.arange(1966,2026)
FIELDS=None


def fields_for_year(year,layer,extra_steps=0):
    arrays=[];times=[];lat=lon=None
    for y,n in [(year,None)]+([(year+1,extra_steps)] if extra_steps else []):
        if not 1966<=y<=2025:raise ValueError('Wind calendar exceeds 1966-2025')
        path=FIELDS/f'{y}.nc'
        if not path.is_file():raise FileNotFoundError(path)
        with xr.open_dataset(path) as ds:
            t=pd.DatetimeIndex(ds.time.values);expected=pd.date_range(f'{y}-01-01',f'{y}-12-31 18:00',freq='6h')
            if not t.equals(expected):raise ValueError('Incomplete six-hourly wind calendar')
            la=ds.latitude.values.astype(float);lo=ds.longitude.values.astype(float)
            if not np.allclose(np.diff(la),-1) or not np.allclose(np.diff(lo),1):raise ValueError('Expected descending 1-degree latitude and ascending longitude')
            if lat is not None and (not np.array_equal(lat,la) or not np.array_equal(lon,lo)):raise ValueError('Wind grids differ between years')
            lat,lon=la,lo
            if n is not None and n>len(t):raise ValueError('Next-year wind prefix exceeds available year')
            sl=slice(None) if n is None else slice(0,n)
            uv=[np.asarray(ds[f'{v}_{layer}'].transpose('time','latitude','longitude').isel(time=sl).values,dtype=np.float32) for v in ('u','v')]
            if not all(np.isfinite(a).all() for a in uv):raise ValueError('Nonfinite wind input')
            arrays.append(uv);times.append(t[sl])
    return np.concatenate([a[0] for a in arrays]),np.concatenate([a[1] for a in arrays]),times[0].append(times[1:]),lat,lon


def prepare_tracks(path):
    t=pd.read_csv(path,parse_dates=['iso_time']);needed={'sid','season','iso_time','lat','lon','wind','nature'}
    if not needed.issubset(t):raise ValueError(f'Missing track columns {needed-set(t)}')
    if t[['sid','season','iso_time','lat','lon','wind']].isna().any().any():raise ValueError('Missing required track values')
    if t.duplicated(['sid','iso_time']).any():raise ValueError('Duplicate SID/time')
    if not t.season.between(1966,2025).all():raise ValueError('Source seasons outside study')
    if not t.iso_time.dt.year.between(1966,2025).all():raise ValueError('Actual timestamps outside ERA5 calendar')
    if not ((t.iso_time.dt.hour%6==0)&(t.iso_time.dt.minute==0)&(t.iso_time.dt.second==0)&(t.iso_time.dt.microsecond==0)&(t.iso_time.dt.nanosecond==0)).all():raise ValueError('Reports must be exactly synoptic')
    if not (t.lat.between(0,40)&t.lon.between(100,180)&t.wind.ge(34)&(t.nature.isna()|t.nature.eq('TS'))).all():raise ValueError('Canonical domain/intensity/nature mismatch')
    if t.groupby('sid').season.nunique().gt(1).any():raise ValueError('Inconsistent source season')
    return t.sort_values(['sid','iso_time']).reset_index(drop=True)

def make_seeds(tracks):
    g=tracks.sort_values(['sid','iso_time']).groupby('sid',sort=True);s=g.first().reset_index()[['sid','season','lat','lon','iso_time']]
    s['report_steps']=g.size().reindex(s.sid).to_numpy();elapsed=(g.iso_time.last()-g.iso_time.first())/pd.Timedelta(hours=6)
    if (elapsed<0).any() or not np.equal(elapsed,np.floor(elapsed)).all() or (elapsed>1464).any():raise ValueError('Invalid elapsed lifetime')
    s['elapsed_steps']=elapsed.reindex(s.sid).to_numpy(int)+1;s['start_year']=s.iso_time.dt.year
    return s

def annual_fields(points):return repo_annual_fields(points,YEARS)

@njit(cache=False)
def integrate(seeds,u,v,lat0,lon0,dy,dx,par,months,mode,nrep=1,capture=True):
    # cols: lat,lon,steps,start_index,source_year_index,realization.
    hist=np.zeros((nrep,60,512));counts=np.zeros(len(seeds),np.int64);bad=0
    maxrecords=int(seeds[:,2].sum()) if capture else 0;paths=np.empty((maxrecords,4));used=0
    for k in range(len(seeds)):
        la,lo=seeds[k,0],seeds[k,1];steps=int(seeds[k,2]);start=int(seeds[k,3]);yy=int(seeds[k,4]);rep=int(seeds[k,5])
        for step in range(steps):
            if not(0<=la<=40 and 100<=lo<=180):break
            j=min(int(la/2.5),15);i=min(int((lo-100)/2.5),31);hist[rep,yy,j*32+i]+=1;counts[k]+=1
            if capture:paths[used]=np.array([k,step,la,lo]);used+=1
            if step==steps-1:break
            if mode==0:idx=start+step
            else:
                mi=start if mode==1 else start+step
                if mi<0 or mi>=len(months):bad+=1;break
                idx=int(months[mi])-1
            if idx<0 or idx>=len(u):bad+=1;break
            su=interp(u[idx],lat0,lon0,dy,dx,la,lo);sv=interp(v[idx],lat0,lon0,dy,dx,la,lo)
            if not(np.isfinite(su) and np.isfinite(sv)):bad+=1;break
            lo+=(par[0]*su+par[1])*21600/(6371000*np.cos(np.deg2rad(la)))*180/np.pi
            la+=(par[0]*sv+par[2])*21600/6371000*180/np.pi
    return hist,paths[:used],counts,bad

def seed_matrix(s,starts,steps,reps=None):return np.column_stack([s.lat,s.lon,steps,starts,s.season.to_numpy()-1966,np.zeros(len(s)) if reps is None else reps]).astype(float)

def path_frame(raw,seeds):
    k=raw[:,0].astype(int);step=raw[:,1].astype(int)
    return pd.DataFrame({'sid':seeds.sid.to_numpy()[k],'season':seeds.season.to_numpy()[k],'iso_time':seeds.iso_time.to_numpy()[k]+pd.to_timedelta(step*6,unit='h'),'lat':raw[:,2],'lon':raw[:,3],'step':step})

def diagnostic_rows(case,layer,annual,paths,obs):
    rows=[]
    for period,sl,ys in [('1966-1995',slice(0,30),(1966,1995)),('1996-2025',slice(30,60),(1996,2025))]:
        field=annual[sl].mean(0);q=paths.loc[paths.season.between(*ys)]
        rows.append({'case':case,'layer':layer,'period':period,'TV_error_vs_observed':float(.5*np.abs(field-obs[sl].mean(0)).sum()),'share_north20':float(field.reshape(16,32)[8:].sum()),'path_points':len(q),**storm_metrics(q)})
    return rows


def simulate_baseline(seeds,layer,par):
    pieces=[]
    for year,seeds_year in seeds.groupby(seeds.iso_time.dt.year,sort=True):
        q=seeds_year.reset_index(drop=True);steps=q.elapsed_steps.to_numpy(int)
        start=pd.Timestamp(f'{year}-01-01');starts=((q.iso_time-start)/pd.Timedelta(hours=6)).to_numpy(int)
        ny=len(pd.date_range(start,f'{year}-12-31 18:00',freq='6h'));extra=max(0,int(np.max(starts+steps))-ny)
        u,v,t,lat,lon=fields_for_year(int(year),layer,extra)
        _,raw,_,bad=integrate(seed_matrix(q,starts,steps),u,v,lat[0],lon[0],lat[1]-lat[0],lon[1]-lon[0],par,t.month.to_numpy(),0,1,True)
        if bad:raise ValueError('Wind input does not cover the complete trajectory')
        pieces.append(path_frame(raw,q))
    paths=pd.concat(pieces,ignore_index=True)
    return annual_fields(paths),paths


def calibrate(tracks,layer):
    pieces=[]
    for _,g in tracks.groupby('sid',sort=False):
        g=g.sort_values('iso_time');la=g.lat.to_numpy();lo=g.lon.to_numpy();t=g.iso_time.to_numpy();dt=np.diff(t)/np.timedelta64(1,'s');keep=(dt>0)&(dt<=43200)
        if not keep.any():continue
        du=np.diff(lo)*np.pi/180*6371000*np.cos(np.deg2rad((la[:-1]+la[1:])/2))/dt;dv=np.diff(la)*np.pi/180*6371000/dt
        q=g.iloc[:-1].loc[keep,['iso_time','lat','lon']].copy();q['du']=du[keep];q['dv']=dv[keep];pieces.append(q)
    q=pd.concat(pieces,ignore_index=True);q['su']=np.nan;q['sv']=np.nan
    for year,g in q.groupby(q.iso_time.dt.year):
        u,v,t,lat,lon=fields_for_year(int(year),layer)
        for r in g.itertuples():
            k=int((r.iso_time-t[0])/pd.Timedelta(hours=6));q.loc[r.Index,['su','sv']]=[interp(f[k],lat[0],lon[0],lat[1]-lat[0],lon[1]-lon[0],r.lat,r.lon) for f in (u,v)]
    data=q[['du','dv','su','sv']].to_numpy()
    if not np.isfinite(data).all():raise ValueError('Nonfinite calibration segment')
    parameters=fit(data)
    return parameters,skill(data,parameters)


def main():
    global FIELDS
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=PROJECT);p.add_argument('--fields',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--refit',action='store_true');a=p.parse_args()
    FIELDS=a.fields.resolve()
    if not all((FIELDS/f'{y}.nc').is_file() for y in YEARS):raise FileNotFoundError('Restore all annual ERA5 inputs described in docs/SUBMONTHLY_STEERING.md')
    a.output.mkdir(parents=True,exist_ok=False)
    tracks=prepare_tracks(a.root/'data/processed/submonthly_steering/domain_ts_standard.csv');seeds=make_seeds(tracks);obs=annual_fields(tracks)
    cal=json.loads((a.root/'results/models/six_hourly/calibration.json').read_text())
    rows=diagnostic_rows('observed','not_applicable',obs,tracks,obs);fields={'observed_annual':obs}
    for layer in ['deep','shallow']:
        if a.refit:
            par,skill_values=calibrate(tracks,layer);cal[layer]['actual_time_parameters']=par.tolist();cal[layer]['actual_skill']=skill_values
        else:par=np.asarray(cal[layer]['actual_time_parameters'])
        annual,paths=simulate_baseline(seeds,layer,par);fields[layer+'_six_hour_actual']=annual
        rows+=diagnostic_rows('six_hour_actual',layer,annual,paths,obs);paths.to_csv(a.output/(layer+'_paths.csv'),index=False)
    pd.DataFrame(rows).to_csv(a.output/'baseline_diagnostics.csv',index=False);np.savez_compressed(a.output/'annual_path_fields.npz',**fields)
    (a.output/'calibration.json').write_text(json.dumps(cal,indent=2))

if __name__=='__main__':main()
