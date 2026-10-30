'Finite beta-advection-model sensitivities and paired year-block uncertainty.'

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import numpy as np

import pandas as pd

import netCDF4 as nc

from numba import njit

from core.paths import PROJECT,output_path,write_json,import_source



@njit(cache=False)

def interp(f,lat0,lon0,dy,dx,lat,lon):

    fy=(lat-lat0)/dy;fx=(lon-lon0)/dx;j=int(np.floor(fy));i=int(np.floor(fx))

    if j<0 or i<0 or j>=f.shape[0]-1 or i>=f.shape[1]-1:return np.nan

    y=fy-j;x=fx-i

    return f[j,i]*(1-y)*(1-x)+f[j,i+1]*(1-y)*x+f[j+1,i]*y*(1-x)+f[j+1,i+1]*y*x



@njit(cache=False)

def simulate(seeds,months,u,v,lat0,lon0,dy,dx,alpha,ub,vb,beta_cos=False):

    hist=np.zeros((60,16,32));points=np.zeros(60,dtype=np.int64)

    for k in range(len(seeds)):

        year=int(seeds[k,0])-1966;lat=seeds[k,1];lon=seeds[k,2];steps=int(seeds[k,3])

        for step in range(steps):

            if not(100<=lon<=180 and 0<=lat<=40):break

            j=min(int(lat/2.5),15);i=min(int((lon-100)/2.5),31)

            hist[year,j,i]+=1;points[year]+=1

            m=months[k,step]-1

            su=interp(u[m],lat0,lon0,dy,dx,lat,lon);sv=interp(v[m],lat0,lon0,dy,dx,lat,lon)

            if not(np.isfinite(su) and np.isfinite(sv)):break

            factor=np.cos(np.deg2rad(lat)) if beta_cos else 1.0

            lon+=(alpha*su+ub*factor)*21600/(6371000*np.cos(np.deg2rad(lat)))*180/np.pi

            lat+=(alpha*sv+vb*factor)*21600/6371000*180/np.pi

    for y in range(60):

        if points[y]>0:hist[y]/=points[y]

    return hist.reshape(60,512),points



def load_fields(path):

    with nc.Dataset(path) as ds:

        tv=ds['time'];dates=nc.num2date(tv[:],tv.units,getattr(tv,'calendar','standard'))

        year=np.array([d.year for d in dates]);month=np.array([d.month for d in dates])

        lat=np.asarray(ds['latitude'][:],float);lon=np.asarray(ds['longitude'][:],float)

        out={}

        for c in ['u','v']:

            keep=np.flatnonzero((year>=1966)&(year<=2025))

            raw=np.asarray(ds[c+'_steer'][keep.min():keep.max()+1],dtype=np.float32)

            yy=year[keep.min():keep.max()+1];mm=month[keep.min():keep.max()+1]

            out[c]=np.stack([raw[(yy>=a)&(yy<=b)&(mm==m)].mean(0,dtype=np.float64)

                            for a,b in [(1966,1995),(1996,2025),(1966,2025)] for m in range(1,13)]).reshape(3,12,len(lat),len(lon))

    return out,lat,lon



def segments(tracks,u,v,lat,lon,period=None,legacy=False):

    lists=[]

    for sid,g in tracks.groupby('sid',sort=False):

        if period and not period[0]<=g.season.iloc[0]<=period[1]:continue

        la=g.lat.to_numpy();lo=g.lon.to_numpy();t=g.iso_time.to_numpy()

        dt=np.diff(t)/np.timedelta64(1,'s')

        ok=(dt>0)&(dt<=43200) if not legacy else dt>0

        divisor=np.full(len(dt),21600.) if legacy else dt

        du=np.diff(lo)*np.pi/180*6371000*np.cos((la[1:]+la[:-1])/2*np.pi/180)/divisor

        dv=np.diff(la)*np.pi/180*6371000/divisor

        su=[];sv=[]

        for j,m in enumerate(g.iso_time.dt.month.to_numpy()[:-1]):

            su.append(interp(u[m-1],lat[0],lon[0],lat[1]-lat[0],lon[1]-lon[0],la[j],lo[j]))

            sv.append(interp(v[m-1],lat[0],lon[0],lat[1]-lat[0],lon[1]-lon[0],la[j],lo[j]))

        su=np.asarray(su);sv=np.asarray(sv)

        ok &=np.isfinite(su)&np.isfinite(sv)&np.isfinite(du)&np.isfinite(dv)

        lists.append(np.column_stack([du[ok],dv[ok],su[ok],sv[ok],np.cos(np.deg2rad(la[:-1][ok]))]))

    return np.concatenate(lists)



def fit(data,zero=False):

    du,dv,su,sv=data[:,:4].T;n=len(data)

    if zero:par=np.array([(su@du+sv@dv)/(su@su+sv@sv),0.,0.])

    else:

        A=np.zeros((n*2,3));A[:n,0]=su;A[n:,0]=sv;A[:n,1]=1;A[n:,2]=1

        par=np.linalg.lstsq(A,np.r_[du,dv],rcond=None)[0]

    return par



def skill(data,par):

    du,dv,su,sv=data[:,:4].T

    return {'n_segments':len(data),'rmse_u_ms':float(np.sqrt(np.mean((du-par[0]*su-par[1])**2))),

        'rmse_v_ms':float(np.sqrt(np.mean((dv-par[0]*sv-par[2])**2)))}



def month_array(seeds,steps,advance):

    maxstep=int(steps.max());a=np.zeros((len(seeds),maxstep),dtype=np.int64)

    for k,(_,r) in enumerate(seeds.iterrows()):

        a[k,:int(steps[k])]=pd.date_range(r.iso_time,periods=int(steps[k]),freq='6h').month if advance else r.iso_time.month

    return a



def run_scenarios(tracks,fields,lat,lon,par,advance=False,elapsed=False,beta_cos=False):

    grouped=tracks.groupby('sid',sort=True)

    seeds=grouped.first().reset_index()[['sid','season','lat','lon','iso_time']]

    count=grouped.size().reindex(seeds.sid).to_numpy()

    if elapsed:

        span=(grouped.iso_time.max()-grouped.iso_time.min()).reindex(seeds.sid)/pd.Timedelta(hours=6)

        count=np.rint(span.to_numpy()).astype(int)+1

    seeds['steps']=count

    e=seeds.loc[seeds.season.le(1995)].copy();l=seeds.loc[seeds.season.ge(1996)].copy()

    def remap(target,donor):

        a=target.copy();q=(a.steps.rank(method='first').to_numpy()-.5)/len(a)

        a['steps']=np.quantile(np.sort(donor.steps.to_numpy()),q).round().astype(int).clip(min=1)

        return a

    scenarios={'early':(e,0),'late':(l,1),'steering':(e,1),'genesis':(remap(l,e),0),'lifetime':(remap(e,l),0)}

    out={}

    for label,(s,p) in scenarios.items():

        steps=s.steps.to_numpy(int);m=month_array(s,steps,advance)

        arr=s[['season','lat','lon','steps']].to_numpy(float)

        h,n=simulate(arr,m,fields['u'][p],fields['v'][p],lat[0],lon[0],lat[1]-lat[0],lon[1]-lon[0],*par,beta_cos)

        select=np.arange(0,30) if label in ['early','steering','lifetime'] else np.arange(30,60)

        if (n[select]==0).any():raise ValueError('A simulation year has no points')

        out[label]=h[select]

    return out



def main():

    bam=import_source('bam_observation_reader','analysis/models/counterfactual_steering.py')

    tracks,_=bam.load_tracks(bam.load_config());years=np.arange(1966,2026)

    xe,ye=bam.grid_edges(bam.load_config(),bam.load_config()['regions']['tc'])

    obs=bam.annual_relative_field(tracks,years,xe,ye).reshape(60,-1)

    ref=obs[30:].mean(0)-obs[:30].mean(0)

    stats=[];baseline=[];calibration=[];allannual={}

    cases=[('850-300',PROJECT/'data/interim/steering.nc'),('850-500',output_path('data/derived/steering_850_500.nc'))]

    for layer,path in cases:

        fields,lat,lon=load_fields(path)

        allseg=segments(tracks,fields['u'][2],fields['v'][2],lat,lon)

        pe=fit(segments(tracks,fields['u'][0],fields['v'][0],lat,lon,(1966,1995)))

        pl=fit(segments(tracks,fields['u'][1],fields['v'][1],lat,lon,(1996,2025)))

        du,dv,su,sv,coslat=allseg.T

        fixed_alpha=(su@(du+coslat)+sv@(dv-2.5*coslat))/(su@su+sv@sv)

        variants=[('actual_dt_max12h',fit(allseg),False,False,False),('zero_intercept',fit(allseg,True),False,False,False),

            ('literature_beta_coslat',np.array([fixed_alpha,-1.,2.5]),False,False,True)]

        if layer=='850-300':variants += [('calibrate_early',pe,False,False,False),('calibrate_late',pl,False,False,False),

            ('calendar_month_updates',fit(allseg),True,False,False),('elapsed_lifetime',fit(allseg),False,True,False)]

        if layer=='850-300':

            import json

            oldpar=json.loads((_LAYOUT_ROOT / 'results/models/counterfactual/criteria.json').read_text())['calibration']

            legacy=run_scenarios(tracks,fields,lat,lon,np.array([oldpar['alpha'],oldpar['u_beta'],oldpar['v_beta']]))

            archived=np.load(_LAYOUT_ROOT / 'results/models/counterfactual/fields.npz')

            mapping={'early':'f_BASE_P1','late':'f_BASE_P2','steering':'f_STEER','genesis':'f_LOC','lifetime':'f_LIFE'}

            err=max(float(abs(legacy[k].mean(0)-archived[v].ravel()).max()) for k,v in mapping.items())

            if err>1e-10:raise AssertionError(f'New engine does not reproduce legacy trajectories: {err}')

            write_json('results/diagnostics/bam_engine_reference_equivalence.json',{'status':'PASS','maxabs_period_fields':err})

        for train,p,validationperiod,fieldperiod in [('early',pe,(1996,2025),1),('late',pl,(1966,1995),0)]:

            val=segments(tracks,fields['u'][fieldperiod],fields['v'][fieldperiod],lat,lon,validationperiod)

            calibration.append({'layer':layer,'trained_period':train,'alpha':p[0],'u_intercept_ms':p[1],'v_intercept_ms':p[2],

                'validation_start':validationperiod[0],'validation_end':validationperiod[1],**skill(val,p)})

        for variant,p,advance,elapsed,beta_cos in variants:

            label=layer+'_'+variant;out=run_scenarios(tracks,fields,lat,lon,p,advance,elapsed,beta_cos)

            for k,v in out.items():allannual[label+'_'+k]=v

            for period,sim,observed in [('early',out['early'].mean(0),obs[:30].mean(0)),('late',out['late'].mean(0),obs[30:].mean(0))]:

                baseline.append({'case':label,'period':period,'r':float(np.corrcoef(sim,observed)[0,1]),'tv_error':float(abs(sim-observed).sum()/2)})

            for scenario,x in [('full',out['late'].mean(0)-out['early'].mean(0)),('steering',out['steering'].mean(0)-out['early'].mean(0)),

                ('genesis',out['genesis'].mean(0)-out['early'].mean(0)),('lifetime',out['lifetime'].mean(0)-out['early'].mean(0))]:

                stats.append({'case':label,'scenario':scenario,'alpha':p[0],'u_intercept_ms':p[1],'v_intercept_ms':p[2],

                    'spatial_r':float(np.corrcoef(x,ref)[0,1]),'tv_amplitude_ratio':float(abs(x).sum()/abs(ref).sum()),

                    'projection_on_observed_axis':float(x@ref/(ref@ref))})

            print('BAM sensitivity completed:',label,flush=True)

    pd.DataFrame(stats).to_csv(_LAYOUT_ROOT / 'results/models/bam_sensitivity_scenarios.csv',index=False)

    pd.DataFrame(baseline).to_csv(_LAYOUT_ROOT / 'results/models/bam_sensitivity_baseline.csv',index=False)

    pd.DataFrame(calibration).to_csv(_LAYOUT_ROOT / 'results/models/bam_cross_period_validation.csv',index=False)

    np.savez(_LAYOUT_ROOT / 'results/models/bam_sensitivity_annual.npz',observed=obs,**allannual)

    # Paired resampling of observed and synthetic years. Intervals are conditional on fitted parameters.

    rng=np.random.default_rng(202406);boot=[]

    def block_draw():

        starts=rng.integers(0,30,10)

        return ((starts[:,None]+np.arange(3))%30).ravel()

    for label in sorted({r['case'] for r in stats}):

        out={k:allannual[label+'_'+k] for k in ['early','late','steering','genesis','lifetime']}

        draws={k:[] for k in ['full','steering','genesis','lifetime']}

        for _ in range(1999):

            i,j=block_draw(),block_draw();d=obs[30:][j].mean(0)-obs[:30][i].mean(0)

            for name,field in [('full',out['late'][j].mean(0)),('steering',out['steering'][i].mean(0)),('genesis',out['genesis'][j].mean(0)),('lifetime',out['lifetime'][i].mean(0))]:

                delta=field-out['early'][i].mean(0);draws[name].append(float(delta@d/(d@d)))

        for name,d in draws.items():

            lo,med,hi=np.quantile(d,[.025,.5,.975]);boot.append({'case':label,'scenario':name,'projection_q025':lo,'projection_median':med,'projection_q975':hi})

    pd.DataFrame(boot).to_csv(_LAYOUT_ROOT / 'results/models/bam_projection_bootstrap.csv',index=False)

    write_json('results/diagnostics/bam_sensitivity.json',{'status':'PASS','reference_reproduction':'separate, exact; historical dt assumption retained only there',

        'new_calibration':'actual segment duration, 0<gap<=12h; joint intercepts are empirical motion residuals, not isolated physical beta drift',

        'uncertainty':'1999 paired circular-three-year bootstrap draws within each period; conditional on calibrated parameters, not full parameter uncertainty',

        'zero_intercept':'diagnostic control, not asserted physically realistic model',

        'literature_beta':'(-1, +2.5) m/s times cos(latitude), adapted parameter-form sensitivity from Lin et al. 2023, doi:10.1029/2023MS003686 Eq.2 and Table 1. Original depth-averaged steering retained; not a full reproduction of their intensity-dependent 850/250-hPa model.',

        'limits':['Scenario fields are not additive causal contributions.','No synoptic winds or storm structure included.',

                  'Finite grouped tests; no optimized parameter sweep or new success threshold.']})



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
