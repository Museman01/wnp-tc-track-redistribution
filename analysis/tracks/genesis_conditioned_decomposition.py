'Genesis-distribution and conditional-track-field decomposition from track records.'

from __future__ import annotations

from core.paths import create_new_output_dir

import numpy as np,pandas as pd



import analysis.shared.record_composition as release_common

import analysis.tracks.genesis_shapley as decomp




def run(catalog, output_dir):

    out=create_new_output_dir(output_dir)

    if catalog=='HKO':

        # HKO native TS reports (IBTrACS 1966-2024 plus official 2025 post-analysis), same domain filter.

        from analysis.shared.agency_catalogues import all_catalogues

        h=all_catalogues()['HKO']['ts_points'];lo,hi,la0,la1=release_common.REGION[0],release_common.REGION[1],release_common.REGION[2],release_common.REGION[3]

        tracks=h.loc[h.season.between(1966,2025)&h.lon.between(lo,hi)&h.lat.between(la0,la1)].sort_values(['season','sid','iso_time'])

    else:tracks=release_common.load_primary_tracks() if catalog=='PRIMARY' else release_common.load_agency_tracks()[catalog]

    summary,boot,fields=decomp.run_one(catalog,tracks,bootstrap=True)

    b=pd.DataFrame(boot)

    for v in ('total_tv','genesis_projection_fraction','propagation_projection_fraction'):

        lo,hi=np.quantile(b[v].dropna(),[0.025,0.975]);summary[f'{v}_boot_lo']=lo;summary[f'{v}_boot_hi']=hi


    np.savez_compressed(out/f'{catalog}_fields.npz',**fields)

    pd.DataFrame([summary]).to_csv(out/f'{catalog}_summary.csv',index=False)

    b.to_csv(out/f'{catalog}_bootstrap.csv',index=False)

    print(catalog,{k:summary[k] for k in summary if 'fraction' in k})



if __name__=='__main__':
    import argparse
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalog',choices=['PRIMARY','USA','TOKYO','CMA','HKO'])
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    if not args.execute:parser.error('add --execute to calculate the decomposition')
    run(args.catalog,args.output_dir)
