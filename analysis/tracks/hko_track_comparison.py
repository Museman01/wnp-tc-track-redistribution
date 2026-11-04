"""Equal-window native four-agency track-field comparison, including common storms."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

from itertools import combinations

import numpy as np

import pandas as pd

from core.paths import PROJECT,output_path,write_json

from analysis.shared.agency_catalogues import all_catalogues, annual_fields, domain

from analysis.shared.inference_stats import compositional_change_test, bh_fdr



def main():

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

    frame.to_csv(_LAYOUT_ROOT / 'results/tracks/hko_track_field_tests.csv',index=False)

    pd.DataFrame(corr).to_csv(_LAYOUT_ROOT / 'results/tracks/hko_difference_correlations.csv',index=False)

    pd.DataFrame({'sid':sorted(common)}).to_csv(_LAYOUT_ROOT / 'results/tracks/four_agency_common_membership.csv',index=False)

    np.savez(_LAYOUT_ROOT / 'results/tracks/hko_track_fields.npz',**arrays)

    write_json('results/diagnostics/hko_comparison.json',{'status':'PASS','window':'1966-2025 (30 early years, 30 late years)',

        'HKO_2025':'Official HKO post-analysis records for 2025, matched to IBTrACS storm identifiers','common_storms':len(common),

        'HKO_flags':'IFLAG character 4; O/V native intensity and O/I/V native position',

        'HKO_native_TS_categories':['TS','STS','T','ST','SuperT'],

        'claim_scope':'Sensitivity to record conventions; shared observations do not imply independent confirmations.',

        'genesis_inference':'Spatial-field agreement alone does not establish small agency or temporal differences in genesis.'})



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
