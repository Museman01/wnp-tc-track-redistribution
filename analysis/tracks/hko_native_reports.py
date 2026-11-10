'Acquire official HKO 2025 post-analysis reports and bind native storms to IBTrACS SIDs.'

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

from urllib.request import urlopen

from datetime import datetime,timezone


import pandas as pd

import numpy as np

from core.paths import PROJECT,output_path,sha256,write_json



def main():

    url='https://data.weather.gov.hk/weatherAPI/hko_data/tc/HKO2025BST.csv'

    target=output_path('data/raw/hko/HKO2025BST.csv')

    if not target.exists():target.write_bytes(urlopen(url,timeout=45).read())

    h=pd.read_csv(target,skiprows=3)

    assert len(h.columns)==12

    h.columns=['name','year','month','day','hour','category','lat100','lon100','pressure','wind','jma_code','hko_code']

    assert h.year.eq(2025).all()

    h['iso_time']=pd.to_datetime(h[['year','month','day']])+pd.to_timedelta(h.hour,unit='h')

    h['lat']=h.lat100/100.;h['lon']=h.lon100/100.

    assert not h.duplicated(['hko_code','iso_time']).any()

    ts_categories=['TS','STS','T','ST','SuperT']

    eligible=set(h.loc[h.category.isin(ts_categories),'hko_code'])

    x=pd.read_csv(PROJECT/'data/raw/ibtracs_wp.csv',skiprows=[1],usecols=['SID','SEASON','NAME','ISO_TIME','LAT','LON'],low_memory=False)

    x=x.loc[x.SEASON.eq(2025)].copy();x['ISO_TIME']=pd.to_datetime(x.ISO_TIME)

    for c in ['LAT','LON']:x[c]=pd.to_numeric(x[c],errors='coerce')

    mapping=[];frames=[]

    for code,g in h.loc[h.hko_code.isin(eligible)].groupby('hko_code'):

        name=g.name.iloc[0].upper().strip();matches=x.loc[x.NAME.str.strip().str.upper().eq(name),'SID'].unique()

        # HKO code 0003 is unnamed; 7 simultaneous reports place it within

        # 39.5 km of this UNNAMED SID, versus >3500 km mean for the only other storm.

        if name=='NAMELESS' and int(code)==3:matches=['2025175N16117']

        if len(matches)!=1:raise ValueError(f'Non-unique HKO/SID binding: {code},{name},{matches}')

        sid=matches[0];g=g.copy();g['sid']=sid;g['season']=2025

        common=g.merge(x.loc[x.SID.eq(sid)],left_on='iso_time',right_on='ISO_TIME')

        if common.empty:raise ValueError('No simultaneous reports for SID verification')

        dy=common.lat-common.LAT;dx=((common.lon-common.LON+180)%360-180)*np.cos(np.deg2rad(common.lat))

        dist=np.sqrt(dx*dx+dy*dy)*111.2

        if dist.max()>500:raise ValueError(f'Investigate position mismatch: {name}, {dist.max()}km')

        mapping.append(dict(hko_code=int(code),name=name,binding_method='explicit unnamed-storm temporal/spatial reconciliation' if name=='NAMELESS' else 'unique same-year storm name plus simultaneous-position verification',jma_code=None if pd.isna(g.jma_code.iloc[0]) else int(g.jma_code.iloc[0]),sid=sid,matched_times=len(common),mean_position_difference_km=float(dist.mean()),max_position_difference_km=float(dist.max())))

        frames.append(g)

    reports=pd.concat(frames,ignore_index=True)

    reports=reports.loc[reports.hour.isin([0,6,12,18])].copy()

    reports['is_native_ts']=reports.category.isin(ts_categories)

    reports['source']='HKO 2025 post-analysis CSV'

    reports.to_csv(output_path('data/derived/hko_2025_native_reports.csv'),index=False)

    pd.DataFrame(mapping).to_csv(_LAYOUT_ROOT / 'results/tracks/hko_2025_sid_binding.csv',index=False)

    receipt=dict(source_url=url,retrieved_utc=datetime.now(timezone.utc).isoformat(),file=str(target.relative_to(output_path('.'))),bytes=target.stat().st_size,sha256=sha256(target),raw_reports=len(h),raw_storms=h.hko_code.nunique(),ts_eligible_storms=len(eligible),standard_time_reports=len(reports),standard_time_ts_reports=int(reports.is_native_ts.sum()),all_eligible_storms_bound_to_unique_frozen_sid=True,coordinate_conversion='native hundredths of degree / 100',native_intensity='HKO category; no wind conversion',history='HKO records for 1966-2024 from IBTrACS and for 2025 from the official post-analysis file')

    write_json('data/requests/hko_2025.json',receipt)

    print('HKO 2025 acquired:',receipt,flush=True)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()

