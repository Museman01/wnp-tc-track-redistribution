"""Equal-window native four-agency track-field comparison, including common storms."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

from itertools import combinations
from pathlib import Path

import numpy as np

import pandas as pd


from analysis.shared.agency_catalogues import all_catalogues, annual_fields, domain

from analysis.shared.inference_stats import compositional_change_test, bh_fdr



def main(output_dir=None):
    out=(_LAYOUT_ROOT / (output_dir or 'results/tracks')).resolve()
    if not out.is_relative_to(_LAYOUT_ROOT.resolve()):raise ValueError('Output directory must remain within the project')
    names=['hko_track_field_tests.csv','hko_difference_correlations.csv','four_agency_common_membership.csv','hko_track_fields.npz']
    if any((out/name).exists() for name in names):raise FileExistsError('Refusing to overwrite existing track results')
    out.mkdir(parents=True,exist_ok=True)
    cats=all_catalogues();years=np.arange(1966,2026)

    common=set.intersection(*(set(domain(c['ts_points']).sid) for c in cats.values()))

    rows=[];corr=[];arrays={'years':years};genesis=[]

    for sample in ['native_all','four_agency_common_storms']:

        ds={}

        for agency,cat in cats.items():

            pts=cat['ts_points'];pts=pts.loc[pts.season.between(1966,2025)]

            if sample!='native_all':pts=pts.loc[pts.sid.isin(common)]

            fields=annual_fields(pts,years)

            r=compositional_change_test(fields,np.arange(30),np.arange(30,60),nperm=9999,block=3,seed=202406)

            ds[agency]=r['change'];arrays[sample+'_'+agency]=r['change'].reshape(16,32)

            rows.append({'sample':sample,'agency':agency,'early_years':30,'late_years':30,

                'n_storms':domain(pts).sid.nunique(),'total_variation':r['tv'],'block_permutation_p':r['global_p'],

                'n_permutations':9999,'block_years':3,'window':'1966-1995 vs 1996-2025'})

            print('HKO comparison',sample,agency,r['tv'],r['global_p'],flush=True)

        for first,second in combinations(ds,2):corr.append({'sample':sample,'agency1':first,'agency2':second,

            'difference_field_r':float(np.corrcoef(ds[first],ds[second])[0,1])})

    frame=pd.DataFrame(rows)

    for sample,idx in frame.groupby('sample').groups.items():frame.loc[idx,'q_bh_four_agency_family']=bh_fdr(frame.loc[idx,'block_permutation_p'])

    frame.to_csv(out / 'hko_track_field_tests.csv',index=False)

    pd.DataFrame(corr).to_csv(out / 'hko_difference_correlations.csv',index=False)

    pd.DataFrame({'sid':sorted(common)}).to_csv(out / 'four_agency_common_membership.csv',index=False)

    np.savez(out / 'hko_track_fields.npz',**arrays)

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    if not args.execute:parser.error('add --execute to calculate track comparisons')
    main(args.output_dir)

