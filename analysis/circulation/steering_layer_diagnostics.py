'Shallow steering, field regressions and spatial circulation diagnostics.'

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import numpy as np

import pandas as pd

import netCDF4 as nc

from scipy.signal import detrend

from core.paths import PROJECT,output_path,write_json

from analysis.shared.agency_catalogues import annual_fields, primary_ts

from analysis.shared.inference_stats import bh_fdr



def coarsen(a):

    ny,nx=(a.shape[-2]//10)*10,(a.shape[-1]//10)*10

    return a[...,:ny,:nx].reshape(*a.shape[:-2],ny//10,10,nx//10,10).mean(axis=(-3,-1))



def field_test(x,y,inference):

    xc=x-x.mean();yc=y-y.mean(axis=0)

    dot=np.tensordot(xc,yc,axes=(0,0));den=np.sqrt(np.sum(xc**2)*np.sum(yc**2,axis=0))

    r=np.divide(dot,den,out=np.zeros_like(dot),where=den>0)

    beta=dot/np.sum(xc**2)

    if not inference:return beta,r,None,None,float(np.mean(r*r)),None

    rng=np.random.default_rng(202406); blocks=[np.arange(i,min(i+3,len(x))) for i in range(0,len(x),3)]

    perms=np.array([np.concatenate([blocks[i] for i in rng.permutation(len(blocks))]) for _ in range(499)])

    dots=xc[perms]@yc.reshape(len(x),-1)

    rr=np.divide(dots,den.ravel()[None,:],out=np.zeros_like(dots),where=den.ravel()[None,:]>0)

    p=((abs(rr)>=abs(r.ravel())[None,:]).sum(0)+1)/500

    pg=(int(np.sum(np.mean(rr*rr,axis=1)>=np.mean(r*r)))+1)/500

    return beta,r,p.reshape(r.shape),bh_fdr(p).reshape(r.shape),float(np.mean(r*r)),pg



def derive(overwrite=False, steering_path=None, output_dir=None):
    if output_dir is None:
        shallow=output_path('data/derived/steering_850_500.nc')
        coarsepath=output_path('data/derived/circulation_monthly_coarse.npz')
    else:
        directory=(PROJECT / output_dir).resolve()
        if not directory.is_relative_to(PROJECT.resolve()):raise ValueError('Output directory must remain within the project')
        directory.mkdir(parents=True,exist_ok=True)
        shallow=directory/'steering_850_500.nc'
        coarsepath=directory/'circulation_monthly_coarse.npz'

    if shallow.exists() and coarsepath.exists():return
    if not overwrite and (shallow.exists() or coarsepath.exists()):
        raise FileExistsError('Refusing to overwrite existing steering inputs')

    with nc.Dataset(PROJECT/'data/interim/era5_wnp_dynamic_plev.nc') as p, nc.Dataset(PROJECT / (steering_path or 'data/interim/steering.nc')) as st:

        lat=np.asarray(p['latitude'][:]);lon=np.asarray(p['longitude'][:]);lev=np.asarray(p['level'][:])

        tv=p['time'];dates=nc.num2date(tv[:],tv.units,getattr(tv,'calendar','standard'))

        keep=np.flatnonzero(np.array([1966<=d.year<=2025 for d in dates]))

        idx=np.array([int(np.where(lev==l)[0][0]) for l in [500,600,700,850]])

        zidx=int(np.where(lev==500)[0][0])

        # Trapezoidal integral over ascending pressure, divided by 350 hPa.

        weights=np.array([50,100,125,75],dtype=float)/350

        coll={k:[] for k in ['u300','v300','u500','v500','z500','u_at500']}

        with nc.Dataset(shallow,'w') as out:

            for dim,n in [('time',len(keep)),('latitude',len(lat)),('longitude',len(lon))]:out.createDimension(dim,n)

            ot=out.createVariable('time','f8',('time',));ot.units=tv.units;ot.calendar=getattr(tv,'calendar','standard');ot[:]=tv[keep]

            out.createVariable('latitude','f4',('latitude',))[:]=lat

            out.createVariable('longitude','f4',('longitude',))[:]=lon

            uv={c:out.createVariable(c+'_steer','f4',('time','latitude','longitude'),zlib=True,complevel=1) for c in ['u','v']}

            for c in uv:uv[c].units='m s-1';uv[c].long_name='pressure-integrated 850-500 hPa steering'

            out.source='ERA5 pressure-level input; derived steering output'

            out.pressure_levels_hpa='500,600,700,850';out.pressure_weights=','.join(map(str,weights))

            for y,year in enumerate(range(1966,2026)):

                ids=keep[y*12:(y+1)*12]

                for c in ['u','v']:

                    raw=np.asarray(p[c][ids,idx,:,:],dtype=np.float64)

                    value=np.einsum('tlji,l->tji',raw,weights).astype('f4')

                    uv[c][y*12:(y+1)*12]=value

                    coll[c+'500'].append(coarsen(value))

                    coll[c+'300'].append(coarsen(np.asarray(st[c+'_steer'][ids],dtype=float)))

                coll['z500'].append(coarsen(np.asarray(p['z'][ids,zidx],dtype=float)/9.80665))

                coll['u_at500'].append(coarsen(np.asarray(p['u'][ids,zidx],dtype=float)))

                if year%10==5:print('Derived shallow steering through',year,flush=True)

    payload={k:np.concatenate(v) for k,v in coll.items()}

    payload.update(year=np.repeat(np.arange(1966,2026),12),month=np.tile(np.arange(1,13),60),

        latitude=lat[:len(lat)//10*10].reshape(-1,10).mean(1),longitude=lon[:len(lon)//10*10].reshape(-1,10).mean(1))

    np.savez(coarsepath,**payload)



def regressions(overwrite=False):
    outputs = [
        'results/tracks/redistribution_index_reconstructed.csv',
        'results/circulation/steering_field_regressions.npz',
        'results/circulation/steering_field_regression_summary.csv',
        'results/circulation/steering_corridor_annual.csv',
        'results/diagnostics/steering_layers.json',
    ]
    if not overwrite and any((_LAYOUT_ROOT / path).exists() for path in outputs):
        raise FileExistsError('Refusing to overwrite existing steering results')

    a=np.load(output_path('data/derived/circulation_monthly_coarse.npz'))

    fields=annual_fields(primary_ts(),np.arange(1966,2026))

    idx=[]

    for t in range(60):

        e=fields[[i for i in range(30) if i!=t]].mean(0)

        l=fields[[i for i in range(30,60) if i!=t]].mean(0)

        d=l-e;m=(l+e)/2;idx.append(2*((fields[t]-m)@d)/(d@d))

    x=np.asarray(idx);x=(x-x.mean())/x.std(ddof=1)

    pd.DataFrame({'year':range(1966,2026),'index_oos_z':x}).to_csv(_LAYOUT_ROOT / 'results/tracks/redistribution_index_reconstructed.csv',index=False)

    seasonal={k:np.array([a[k][(a['year']==year)&np.isin(a['month'],[6,7,8,9,10])].mean(0) for year in range(1966,2026)]) for k in ['u300','v300','u500','v500','z500','u_at500']}

    z=seasonal['z500'];ze=z-z.mean(axis=2,keepdims=True)

    res={'latitude':a['latitude'],'longitude':a['longitude']};rows=[]

    for layer in ['300','500']:

        for scale in ['raw','detrended']:

            xx=x if scale=='raw' else detrend(x)

            for name,y in [('u',seasonal['u'+layer]),('v',seasonal['v'+layer]),('z_eddy',ze)]:

                yy=y if scale=='raw' else detrend(y,axis=0)

                beta,r,p,q,r2,pg=field_test(xx,yy,scale=='detrended')

                res[f'{layer}_{scale}_{name}_beta']=beta

                rows.append({'layer_hpa':f'850-{layer}','timescale':scale,'variable':name,'global_mean_r2':r2,

                    'global_block_permutation_p':pg,'n_cells_fdr_q_lt_0_05':None if q is None else int((q<.05).sum()),'max_abs_r':float(abs(r).max())})

    # Period-difference geometry is descriptive; regression covariance is kept separate.

    for k,value in seasonal.items():res['period_change_'+k]=value[30:].mean(0)-value[:30].mean(0)

    res['period_change_z_eddy']=ze[30:].mean(0)-ze[:30].mean(0)

    res['u500_early']=seasonal['u_at500'][:30].mean(0);res['u500_late']=seasonal['u_at500'][30:].mean(0)

    res['z500_early']=z[:30].mean(0);res['z500_late']=z[30:].mean(0)

    np.savez(_LAYOUT_ROOT / 'results/circulation/steering_field_regressions.npz',**res)

    out=pd.DataFrame(rows);out.to_csv(_LAYOUT_ROOT / 'results/circulation/steering_field_regression_summary.csv',index=False)

    write_json('results/diagnostics/steering_layers.json',{

        'pressure_weights_500_600_700_850':[50/350,100/350,125/350,75/350],

        'inference':'499 ordered-three-year block permutations; detrended fields; local BH FDR',

        'geometry':'descriptive two-period fields, not causal attribution','coverage':'1966-2025; annual June-October; original coarsening retained'})

    # Same unweighted corridor averaging as the original Figure 5c source.

    la=a['latitude'];lo=a['longitude'];lm=(la>=10)&(la<=30);om=(lo>=110)&(lo<=160)

    corridor=[]

    for layer in ['300','500']:

        for yi,year in enumerate(range(1966,2026)):

            corridor.append({'layer_hpa':f'850-{layer}','year':year,'index_oos_z':x[yi],

                'u_ms':float(seasonal['u'+layer][yi][lm][:,om].mean()),

                'v_ms':float(seasonal['v'+layer][yi][lm][:,om].mean())})

    pd.DataFrame(corridor).to_csv(_LAYOUT_ROOT / 'results/circulation/steering_corridor_annual.csv',index=False)




if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    overwrite='--overwrite' in sys.argv
    derive(overwrite=overwrite);regressions(overwrite=overwrite)

