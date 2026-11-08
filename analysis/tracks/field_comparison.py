"""Compare storm-equal shared-genesis-bin and full track-point spatial contrasts."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from core.paths import PROJECT

def compare(root):
    rows=[]
    with np.load(root/'data/reference/environment_validation/supplemental_validation/results/wnp_tc_redistribution_pattern.npz',allow_pickle=True) as full, np.load(root/'results/tracks/hko_track_fields.npz') as hko:
        years=full['years']
        for record in ['PRIMARY','USA','TOKYO','CMA','HKO']:
            if record=='HKO':point=hko['native_all_HKO'].ravel()
            else:
                annual=full[record+'_track_point_annual_fields'];point=annual[years>=1996].mean(0)-annual[years<=1995].mean(0)
            with np.load(root/f'results/genesis/{record}_fields.npz') as fields:
                total=fields['total_field'].ravel();conditional=fields['propagation_field'].ravel()
                rows.append(dict(record=record,spatial_r=float(np.corrcoef(total,point)[0,1]),projection=float(total@point/(point@point)),conditional_share=float(conditional@total/(total@total))))
    return pd.DataFrame(rows)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=PROJECT);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    compare(a.root).to_csv(a.output,index=False)
