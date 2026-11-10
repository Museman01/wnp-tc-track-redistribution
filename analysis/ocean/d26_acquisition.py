"""Submit and collect CDS ORAS5 D26 requests."""

from __future__ import annotations

import argparse

import json

import re

from cdsapi.api import get_url_key_verify

from ecmwf.datastores import Client

from core.paths import output_path, write_json, sha256



def main():

    parser=argparse.ArgumentParser()

    parser.add_argument('action',choices=['submit-series','collect-series'])

    args=parser.parse_args()

    try:

        url,key,verify=get_url_key_verify(None,None,None)

        client=Client(url=url,key=key,verify=verify,timeout=45,maximum_tries=1,sleep_max=5,progress=False,log_callback=lambda *a,**k:None)


        # One variable only, all calendar months, same full comparison as OHC300.

        # Five consolidated batches plus one operational batch; no ensemble or 3-D profiles.

        batches=[(start,min(start+9,2014),'consolidated') for start in range(1966,2015,10)]

        batches += [(2015,2025,'operational')]

        for start,end,product in batches:

            rel=f'data/requests/d26_{start}_{end}.json'; path=output_path(rel)

            if args.action=='submit-series':

                if path.exists():print('Retained existing request',start,end,flush=True);continue

                request={'product_type':[product],'vertical_resolution':'single_level',

                    'variable':['depth_of_26_c_isotherm'],'year':[str(y) for y in range(start,end+1)],

                    'month':[f'{m:02d}' for m in range(1,13)]}

                remote=client.submit('reanalysis-oras5',request)

                write_json(rel,{'collection':'reanalysis-oras5','request_id':remote.request_id,'request':request,'status':remote.status})

                print('Submitted D26',start,end,remote.status,flush=True)

            else:

                if not path.exists():continue

                info=json.loads(path.read_text(encoding='utf-8'))

                target=output_path(f'data/raw/d26/oras5_{product}_d26_{start}_{end}.zip')

                if info.get('download_complete') and target.exists():print('Retained complete download',start,end,flush=True);continue

                remote=client.get_remote(info['request_id']);info['status']=remote.status

                write_json(rel,info)

                print('D26',start,end,info['status'],flush=True)

                if info['status']=='successful':

                    remote.download(str(target))

                    info.update(download_complete=True,file=str(target.relative_to(output_path('.'))),bytes=target.stat().st_size,sha256=sha256(target))

                    write_json(rel,info)

                    print('Downloaded',start,end,target.stat().st_size,flush=True)


    except Exception as exc:

        message=re.sub(r'(?i)(key|token|authorization)([=: ]+)[^\s,;]+',r'\1\2[REDACTED]',str(exc))

        write_json('data/requests/d26_acquisition_error.json',{'error_type':type(exc).__name__,'message':message,'action':args.action})

        print('D26 request unavailable:',type(exc).__name__,message[:1200])

        raise SystemExit(2)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()

