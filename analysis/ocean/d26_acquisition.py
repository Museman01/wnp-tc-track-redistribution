"""Submit or inspect a single bounded CDS ORAS5 D26 probe; never print credentials."""

from __future__ import annotations

import argparse

import json

import re

from cdsapi.api import get_url_key_verify

from ecmwf.datastores import Client

from core.paths import output_path, write_json, sha256



def main():

    parser=argparse.ArgumentParser()

    parser.add_argument('action',choices=['submit-probe','status','download-probe','submit-series','collect-series'])

    args=parser.parse_args()

    receipt=output_path('data/requests/d26_probe.json')

    try:

        url,key,verify=get_url_key_verify(None,None,None)

        client=Client(url=url,key=key,verify=verify,timeout=45,maximum_tries=1,sleep_max=5,progress=False,log_callback=lambda *a,**k:None)

        if args.action in ['submit-series','collect-series']:

            # One variable only, all calendar months, same full comparison as OHC300.

            # Six consolidated batches plus one operational batch; no ensemble or 3-D profiles.

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

            return

        if args.action=='submit-probe':

            if receipt.is_file():

                print('Existing probe receipt retained; use status.'); return

            request={'product_type':['consolidated'],'vertical_resolution':'single_level',

                     'variable':['depth_of_26_c_isotherm'],'year':['2014'],'month':['07']}

            write_json('data/requests/d26_probe_request.json',request)

            remote=client.submit('reanalysis-oras5',request)

            payload={'collection':'reanalysis-oras5','request_id':remote.request_id,'request':request,

                     'status':remote.status,'catalogue':'https://cds.climate.copernicus.eu/datasets/reanalysis-oras5',

                     'scope':'one global single-level month; inspect semantics before full retrieval'}

            write_json('data/requests/d26_probe.json',payload)

            print('D26 probe submitted:',remote.request_id,remote.status)

        else:

            info=json.loads(receipt.read_text(encoding='utf-8'))

            remote=client.get_remote(info['request_id'])

            state=remote.status

            info['status']=state

            write_json('data/requests/d26_probe.json',info)

            print('D26 probe status:',state)

            if args.action=='download-probe' and state=='successful':

                target=output_path('data/raw/d26/oras5_consolidated_d26_201407.zip')

                if not target.exists(): remote.download(str(target))

                write_json('data/requests/d26_probe_download.json',{'file':str(target),'bytes':target.stat().st_size,'sha256':sha256(target)})

                print('Probe downloaded:',target.stat().st_size,'bytes')

    except Exception as exc:

        message=re.sub(r'(?i)(key|token|authorization)([=: ]+)[^\s,;]+',r'\1\2[REDACTED]',str(exc))

        write_json('data/requests/d26_acquisition_error.json',{'error_type':type(exc).__name__,'message':message,'action':args.action})

        print('D26 probe unavailable:',type(exc).__name__,message[:1200])

        raise SystemExit(2)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
