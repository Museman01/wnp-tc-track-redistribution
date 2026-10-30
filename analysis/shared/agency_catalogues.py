'Agency catalogues.'
from __future__ import annotations

import numpy as np

import pandas as pd

from core.paths import PROJECT,output_path

from analysis.shared.agency_data import read_ibtracs_agencies, build_agency_catalog



def primary_ts():

    x=pd.read_csv(PROJECT/'data/processed/tracks.csv',parse_dates=['iso_time'])

    return x.loc[x.season.between(1966,2025)&x.wind.ge(34)&(x.nature.eq('TS')|x.nature.isna())].copy()



def all_catalogues(start=1966):

    frame=read_ibtracs_agencies(PROJECT/'data/raw/ibtracs_wp.csv',start,2025)

    cats={key:build_agency_catalog(frame,key) for key in ['USA','TOKYO','CMA']}

    h=pd.read_csv(PROJECT/'data/raw/ibtracs_wp.csv',skiprows=[1],usecols=[

        'SID','SEASON','ISO_TIME','TRACK_TYPE','IFLAG','HKO_LAT','HKO_LON','HKO_CAT','HKO_WIND'],low_memory=False)

    for c in ['SEASON','HKO_LAT','HKO_LON','HKO_WIND']:h[c]=pd.to_numeric(h[c],errors='coerce')

    h['ISO_TIME']=pd.to_datetime(h.ISO_TIME,errors='coerce')

    h=h.loc[h.SEASON.between(max(start,1961),2024)&h.ISO_TIME.dt.hour.isin([0,6,12,18])&~h.TRACK_TYPE.fillna('').str.lower().str.contains('spur')].copy()

    flag=h.IFLAG.fillna('').str.get(3)

    ts=flag.isin(['O','V'])&h.HKO_CAT.str.strip().isin(['TS','STS','T','ST','SuperT'])&h.HKO_LAT.notna()&h.HKO_LON.notna()

    eligible=set(h.loc[ts,'SID'])

    pos=flag.isin(['O','I','V'])&h.SID.isin(eligible)&h.HKO_LAT.notna()&h.HKO_LON.notna()

    def convert(mask):

        x=h.loc[mask,['SID','SEASON','ISO_TIME','HKO_LAT','HKO_LON','HKO_WIND']].copy()

        x.columns=['sid','season','iso_time','lat','lon','wind']

        x['lon']%=360

        return x.sort_values(['sid','iso_time']).drop_duplicates(['sid','iso_time'],keep='last')

    cats['HKO']={'agency':'HKO','ts_points':convert(ts),'tracks':convert(pos),'eligible_sids':eligible}

    addition=output_path('data/derived/hko_2025_native_reports.csv')

    if addition.exists():

        a=pd.read_csv(addition,parse_dates=['iso_time'])

        columns=['sid','season','iso_time','lat','lon','wind']

        cats['HKO']['ts_points']=pd.concat([cats['HKO']['ts_points'],a.loc[a.is_native_ts,columns]],ignore_index=True)

        cats['HKO']['tracks']=pd.concat([cats['HKO']['tracks'],a[columns]],ignore_index=True)

        cats['HKO']['eligible_sids']|=set(a.sid)

    return cats



def domain(x, north=40):

    return x.loc[x.lon.between(100,180)&x.lat.between(0,north)].copy()



def annual_fields(points,years):

    xe=np.arange(100,182.5,2.5); ye=np.arange(0,42.5,2.5)

    rows=[]

    for year in years:

        a=domain(points.loc[points.season.eq(year)])

        h=np.histogram2d(a.lon,a.lat,bins=[xe,ye])[0].T

        if h.sum()==0:raise ValueError(f'Missing spatial data in {year}; never insert zero')

        rows.append((h/h.sum()).ravel())

    return np.asarray(rows)
