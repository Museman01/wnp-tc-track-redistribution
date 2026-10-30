'Deep- and shallow-layer monthly storm-sampling decomposition.'

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import numpy as np

import pandas as pd

import netCDF4 as nc

from core.paths import PROJECT,output_path,write_json

from analysis.shared.agency_catalogues import all_catalogues



def coarsen_file(path):

    out={c:[] for c in ['u','v']}

    with nc.Dataset(path) as ds:

        lat=np.asarray(ds['latitude'][:]);lon=np.asarray(ds['longitude'][:])

        la=np.where((lat>=0)&(lat<=40))[0];lo=np.where((lon>=100)&(lon<=180))[0]

        tv=ds['time'];dates=nc.num2date(tv[:],tv.units,getattr(tv,'calendar','standard'))

        yy=np.array([d.year for d in dates]);mm=np.array([d.month for d in dates])

        lat=lat[la];lon=lon[lo]

        jl=np.clip(np.digitize(lat,np.arange(0,42.5,2.5))-1,0,15)

        il=np.clip(np.digitize(lon,np.arange(100,182.5,2.5))-1,0,31)

        A=np.zeros((16,len(lat)));A[jl,np.arange(len(lat))]=np.cos(np.deg2rad(lat))

        B=np.zeros((32,len(lon)));B[il,np.arange(len(lon))]=1

        den=A.sum(1)[:,None]*B.sum(1)[None,:]

        for year in range(1966,2026):

            ids=np.flatnonzero((yy==year)&np.isin(mm,[6,7,8,9,10]))

            if len(ids)!=5:raise ValueError('Missing steering months')

            for c in ['u','v']:

                raw=np.asarray(ds[c+'_steer'][ids,la.min():la.max()+1,lo.min():lo.max()+1],dtype=float)

                out[c].append(np.einsum('jJ,tJI,iI->tji',A,raw,B,optimize=True)/den)

    return {c:np.stack(v).reshape(60,-1) for c,v in out.items()}



def distributions(tracks):

    n=np.zeros((60,5,16,32));t=tracks.loc[tracks.season.between(1966,2025)&tracks.iso_time.dt.month.isin([6,7,8,9,10])&tracks.lon.between(100,180)&tracks.lat.between(0,40)]

    y=t.season.to_numpy(int)-1966;m=t.iso_time.dt.month.to_numpy(int)-6

    j=np.clip(np.digitize(t.lat,np.arange(0,42.5,2.5))-1,0,15)

    i=np.clip(np.digitize(t.lon,np.arange(100,182.5,2.5))-1,0,31)

    np.add.at(n,(y,m,j,i),1)

    if (n.sum((1,2,3))==0).any():raise ValueError('No annual storm-sampling data')

    return n/n.sum((1,2,3),keepdims=True)



def values(C,weights1,weights2):

    a=np.einsum('bi,ij,bj->b',weights1,C,weights1,optimize=True)

    b=np.einsum('bi,ij,bj->b',weights1,C,weights2,optimize=True)

    c=np.einsum('bi,ij,bj->b',weights2,C,weights1,optimize=True)

    d=np.einsum('bi,ij,bj->b',weights2,C,weights2,optimize=True)

    delta=(weights2-weights1)@np.diag(C)

    field=.5*(d-c+b-a);sampling=.5*(d+c-b-a);cos=delta-d+a

    return {'total_change':delta,'field':field,'sampling':sampling,'cosampling':cos}



def main():

    rows=[];checks=[];cats=all_catalogues()

    rng=np.random.default_rng(202406);blocks=np.arange(60).reshape(20,3)

    first=np.zeros((10000,60));second=np.zeros_like(first)

    first[0,:30]=1/30;second[0,30:]=1/30

    for k in range(1,10000):

        order=rng.permutation(20)

        first[k,blocks[order[:10]].ravel()]=1/30;second[k,blocks[order[10:]].ravel()]=1/30

    for layer,path in [('850-300',PROJECT/'data/interim/steering.nc'),('850-500',output_path('data/derived/steering_850_500.nc'))]:

        cache=output_path('data/derived/sampling_'+layer+'.npz')

        if cache.exists():U=dict(np.load(cache))

        else:U=coarsen_file(path);np.savez(cache,**U)

        for agency in ['USA','TOKYO','CMA']:

            n=distributions(cats[agency]['tracks']);N=n.reshape(60,-1)

            pool=n.mean(0);sp=pool/pool.sum((1,2),keepdims=True)

            nf=n.sum((2,3))[:,:,None,None]*sp[None,:,:,:]

            for c in ['u','v']:

                if not np.isfinite(U[c]).all():raise ValueError('Nonfinite monthly field requires explicit support handling')

                C=N@U[c].T;F=nf.reshape(60,-1)@U[c].T

                a=values(C,first,second);f=values(F,first,second)

                residual=float(abs(a['total_change']-a['field']-a['sampling']-a['cosampling']).max())

                if residual>1e-10:raise AssertionError(residual)

                checks.append({'layer':layer,'agency':agency,'component':c,'maxabs_additivity':residual})

                for name,v in {**a,'fixed_field':f['field'],'fixed_seasonal':f['sampling'],'fixed_cosampling':f['cosampling']}.items():

                    p=(int((abs(v[1:])>=abs(v[0])).sum())+1)/10000

                    rows.append({'layer':layer,'agency':agency,'wind_component':c,'term':name,'change_ms':float(v[0]),

                        'two_sided_block_permutation_p':p,'n_permutations':9999,'block_years':3})

            print('Steering sampling',layer,agency,'complete',flush=True)

    pd.DataFrame(rows).to_csv(_LAYOUT_ROOT / 'results/circulation/steering_sampling_layers.csv',index=False)

    write_json('results/diagnostics/steering_sampling_layers.json',{'status':'PASS','closure_checks':checks,

        'null':'9999 ordered-three-year block permutations, 30/30 split; paired annual field and occupancy labels',

        'fixed_space':'pooled monthly spatial distribution, invariant to equal-size period relabelling',

        'meaning':'sampling includes the explained track-position change; not a causal mechanism',

        'design':'layer sensitivity added after original results, not prespecified historical analysis'})



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
