'D26 time coordinates.'
from __future__ import annotations

import json

import re

import zipfile

import hashlib

import numpy as np

import pandas as pd

import netCDF4 as nc

from core.paths import output_path,write_json



def main():

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

    checked=[]

    for year in range(1966,2026):

        with nc.Dataset(output_path(f'data/derived/d26/oras5_d26_{year}.nc'),'r+') as d:

            before=hashlib.sha256(np.ma.filled(d['d26'][:],np.nan).tobytes()).hexdigest()

            t=d['time'];old=nc.num2date(t[:],t.units,t.calendar)

            assert [(x.year,x.month) for x in old]==[(year,m) for m in range(1,13)]

            t[:]=nc.date2num([times[year,m] for m in range(1,13)],t.units,t.calendar)

            t.timestamp_provenance='Native ORAS5 time_counter, preserved from each CDS archive member.'

            after=hashlib.sha256(np.ma.filled(d['d26'][:],np.nan).tobytes()).hexdigest()

            assert before==after

            checked.append(dict(year=year,data_sha256=after,changed_times=sum(str(a)!=str(times[year,m]) for m,a in enumerate(old,1))))

    pd.DataFrame(rows).sort_values(['year','month']).to_csv(output_path('data/derived/d26_native_times.csv'),index=False)

    write_json('results/diagnostics/d26_native_times.json',dict(status='PASS',source_months=720,years=checked,numerical_effect='None: D26 values and calendar-month ordering unchanged; native source midpoint dates restored.'))



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
