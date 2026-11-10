'D26 time coordinates.'
from __future__ import annotations

import json

import re

import zipfile

import hashlib
import shutil
from pathlib import Path

import numpy as np

import pandas as pd

import netCDF4 as nc

from core.paths import PROJECT,output_path



def main(output_dir):
    out=(PROJECT / output_dir).resolve()
    if not out.is_relative_to(PROJECT.resolve()):raise ValueError('Output directory must remain within the project')
    if out.exists():raise FileExistsError(f'Output directory already exists: {out}')
    times={};rows=[]

    for receipt in sorted(output_path('data/requests').glob('d26_*_*.json')):

        info=json.loads(receipt.read_text(encoding='utf-8'))

        if not info.get('download_complete'):continue

        archive=output_path(info['file'])

        if re.search(r'_\d{6}\.zip$',archive.name):continue

        with zipfile.ZipFile(archive) as z:

            for name in sorted(z.namelist()):

                if not name.endswith('.nc'):continue

                year,month=map(int,re.search(r'_(\d{4})(\d{2})_',name).groups())

                with nc.Dataset('memory',memory=z.read(name)) as d:

                    t=d['time_counter'];calendar=getattr(t,'calendar','standard')

                    stamp=nc.num2date(t[0],t.units,calendar)

                    assert (stamp.year,stamp.month)==(year,month)

                    assert (year,month) not in times

                    times[year,month]=stamp

                    rows.append(dict(year=year,month=month,source_timestamp=str(stamp),source_units=t.units,source_calendar=calendar,archive=archive.name,member=name))

        print('Source timestamps read:',archive.name,flush=True)

    assert len(times)==720

    out.mkdir(parents=True)
    for year in range(1966,2026):
        source=output_path(f'data/derived/d26/oras5_d26_{year}.nc')
        target=out/source.name
        shutil.copy2(source,target)
        with nc.Dataset(target,'r+') as d:

            before=hashlib.sha256(np.ma.filled(d['d26'][:],np.nan).tobytes()).hexdigest()

            t=d['time'];old=nc.num2date(t[:],t.units,t.calendar)

            assert [(x.year,x.month) for x in old]==[(year,m) for m in range(1,13)]

            t[:]=nc.date2num([times[year,m] for m in range(1,13)],t.units,t.calendar)

            t.timestamp_provenance='Native ORAS5 time_counter, preserved from each CDS archive member.'

            after=hashlib.sha256(np.ma.filled(d['d26'][:],np.nan).tobytes()).hexdigest()

            assert before==after

    pd.DataFrame(rows).sort_values(['year','month']).to_csv(out/'d26_native_times.csv',index=False)



if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    if not args.execute:parser.error('add --execute to prepare updated files')
    main(args.output_dir)

