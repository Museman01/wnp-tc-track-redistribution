'D26 preparation.'
from __future__ import annotations

import zipfile

import re

import json
from pathlib import Path

import numpy as np

import pandas as pd

import netCDF4 as nc

from core.paths import PROJECT,output_path



def main(output_dir):
    output_root=(PROJECT / output_dir).resolve()
    if not output_root.is_relative_to(PROJECT.resolve()):raise ValueError('Output directory must remain within the project')
    if output_root.exists():raise FileExistsError(f'Output directory already exists: {output_root}')
    files=[]

    for receipt in output_path('data/requests').glob('d26_*_*.json'):

        info=json.loads(receipt.read_text(encoding='utf-8'))

        if info.get('download_complete'):files.append(output_path(info['file']))

    files=sorted(files)
    if not files:raise FileNotFoundError('No downloaded D26 archives')
    output_root.mkdir(parents=True)
    inventory=[]; pending={}; pending_times={}; reference=None

    for archive in files:

        if re.search(r'_\d{6}\.zip$',archive.name):continue # skip single-month archives

        with zipfile.ZipFile(archive) as z:

            for name in sorted(z.namelist()):

                if not name.endswith('.nc'):continue

                m=re.search(r'_(\d{4})(\d{2})_',name)

                if not m:raise ValueError(f'No date in ORAS5 member: {name}')

                year,month=map(int,m.groups())

                if not 1966<=year<=2025:raise ValueError(year)

                final=output_root/f'oras5_d26_{year}.nc'

                if final.exists():continue

                with nc.Dataset('memory',memory=z.read(name)) as d:

                    if d['so26chgt'].units!='m':raise ValueError('Unexpected D26 units')

                    native_time=d['time_counter']

                    timestamp=nc.num2date(native_time[0],native_time.units,getattr(native_time,'calendar','standard'))

                    if (timestamp.year,timestamp.month)!=(year,month):raise ValueError('Source month mismatch')

                    lat=np.asarray(d['nav_lat'][:]);lon=np.asarray(d['nav_lon'][:])%360

                    mask=(lat>=0)&(lat<=50)&(lon>=100)&(lon<=180)

                    ij=np.argwhere(mask); ymin,xmin=ij.min(0);ymax,xmax=ij.max(0)+1

                    submask=mask[ymin:ymax,xmin:xmax]

                    la=lat[ymin:ymax,xmin:xmax];lo=lon[ymin:ymax,xmin:xmax]

                    arr=np.ma.filled(d['so26chgt'][0,ymin:ymax,xmin:xmax],np.nan).astype('f4')

                    arr[~submask]=np.nan

                    if np.nanmin(arr)<0:raise ValueError('Negative depth needs investigation')

                    if reference is None:reference=(la.copy(),lo.copy())

                    elif not(np.allclose(la,reference[0],atol=5e-5)&np.allclose(lo,reference[1],atol=5e-5)):raise ValueError('Grid drift')

                    inventory.append({'year':year,'month':month,'archive':archive.name,'member':name,

                        'finite_ocean_cells':int(np.isfinite(arr).sum()),'native_zero_cells':int((arr==0).sum()),

                        'missing_cells':int(np.isnan(arr).sum()),'min_m':float(np.nanmin(arr)),'max_m':float(np.nanmax(arr))})

                    pending.setdefault(year,{})[month]=arr

                    pending_times.setdefault(year,{})[month]=timestamp

                    if len(pending[year])==12:

                        with nc.Dataset(final,'w') as out:

                            for dim,size in [('time',12),('y',len(la)),('x',la.shape[1])]:out.createDimension(dim,size)

                            tv=out.createVariable('time','f8',('time',));tv.units='days since 1960-01-01';tv.calendar='standard'

                            tv[:]=nc.date2num([pending_times[year][k] for k in range(1,13)],tv.units,calendar=tv.calendar)

                            for var,value in [('nav_lat',la),('nav_lon',lo),('region_mask',submask.astype('i1'))]:

                                out.createVariable(var,'i1' if var=='region_mask' else 'f4',('y','x'))[:]=value

                            v=out.createVariable('d26','f4',('time','y','x'),zlib=True,complevel=2,fill_value=np.nan)

                            v.units='m';v.long_name='ORAS5 monthly native depth of 26C isotherm'

                            v[:]=np.stack([pending[year][k] for k in range(1,13)])

                            out.zero_handling='Native delivered zeros preserved; missing values are not converted to zero.'

                            out.production_stream='consolidated' if year<=2014 else 'operational'

                        del pending[year]

                        del pending_times[year]


    frame=pd.DataFrame(inventory)
    if len(frame):frame.sort_values(['year','month']).to_csv(output_root/'d26_month_inventory.csv',index=False)
    present=[y for y in range(1966,2026) if (output_root/f'oras5_d26_{y}.nc').is_file()]
    if len(present)!=60 or pending:
        raise ValueError(f'Incomplete D26 input: {len(present)} of 60 years complete')


if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    if not args.execute:parser.error('add --execute to prepare D26 fields')
    main(args.output_dir)

