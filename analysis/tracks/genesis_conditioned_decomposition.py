'Genesis-distribution and conditional-track-field decomposition from track records.'

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import importlib.util,sys,types,json

import numpy as np,pandas as pd

from core.paths import PROJECT,output_path,write_json



import analysis.shared.record_composition as release_common

import analysis.tracks.genesis_shapley as decomp



OUT=_LAYOUT_ROOT / 'results/genesis'

# Manuscript values (Section 3.3): conditional share and 95% bootstrap interval, genesis share.

EXPECTED={'PRIMARY':(79.1,68.8,93.5,20.9),'USA':(79.4,69.4,93.6,20.6),'TOKYO':(73.7,65.4,92.3,26.3),'CMA':(80.0,67.6,93.6,20.0)}



def run(catalog):

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

    OUT.mkdir(parents=True,exist_ok=True)

    np.savez_compressed(OUT/f'{catalog}_fields.npz',**fields)

    pd.DataFrame([summary]).to_csv(OUT/f'{catalog}_summary.csv',index=False)

    b.to_csv(OUT/f'{catalog}_bootstrap.csv',index=False)

    print(catalog,{k:summary[k] for k in summary if 'fraction' in k})



def check():

    rows=[];ok=True

    for cat,(cond,lo,hi,gen) in EXPECTED.items():

        s=pd.read_csv(OUT/f'{cat}_summary.csv').iloc[0]

        got=(round(100*s.propagation_projection_fraction,1),round(100*s.propagation_projection_fraction_boot_lo,1),

             round(100*s.propagation_projection_fraction_boot_hi,1),round(100*s.genesis_projection_fraction,1))

        match=got==(cond,lo,hi,gen);ok&=match

        rows.append(dict(catalog=cat,expected=[cond,lo,hi,gen],regenerated=list(got),match=match))

    write_json('results/diagnostics/genesis_decomposition_regeneration.json',dict(status='PASS' if ok else 'MISMATCH',rows=rows,

        code='analysis/tracks/genesis_shapley.py with the supplied input records'))

    print('PASS' if ok else 'MISMATCH',rows)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    check() if sys.argv[1]=='check' else run(sys.argv[1])
