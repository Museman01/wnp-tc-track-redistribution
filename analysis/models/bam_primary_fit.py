'Bam primary fit.'
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import sys

import numpy as np, pandas as pd

from core.paths import PROJECT, PROJECT_ROOT, output_path, import_source



def main():

    import analysis.models.bam_sensitivity as bs

    fields, lat, lon = bs.load_fields(PROJECT/'data/interim/steering.nc')

    bam = import_source('bam_observation_reader', 'analysis/models/counterfactual_steering.py')

    tracks, _ = bam.load_tracks(bam.load_config())

    seg = bs.segments(tracks, fields['u'][2], fields['v'][2], lat, lon)

    legacy = bs.segments(tracks, fields['u'][2], fields['v'][2], lat, lon, legacy=True)

    rows = []

    for name, data in [('elapsed_time_primary', seg), ('fixed_6h_steps', legacy)]:

        par = bs.fit(data); sk = bs.skill(data, par)

        du, dv, su, sv = data[:, :4].T

        rows.append(dict(calibration=name, n_pairs=len(data), alpha=par[0], u_intercept_ms=par[1], v_intercept_ms=par[2],

                         r_u=float(np.corrcoef(du, par[0]*su+par[1])[0, 1]), r_v=float(np.corrcoef(dv, par[0]*sv+par[2])[0, 1]),

                         rmse_u_ms=sk['rmse_u_ms'], rmse_v_ms=sk['rmse_v_ms']))

    # summed single-factor residual for the primary calibration

    a = np.load(_LAYOUT_ROOT / 'results/models/bam_sensitivity_annual.npz')

    m = {k: a['850-300_actual_dt_max12h_'+k].mean(0) for k in ['early', 'late', 'steering', 'genesis', 'lifetime']}

    d = {k: m[k]-m['early'] for k in ['late', 'steering', 'genesis', 'lifetime']}

    resid = d['late']-(d['steering']+d['genesis']+d['lifetime'])

    old = np.load(_LAYOUT_ROOT / 'results/models/counterfactual/fields.npz')

    resid_old = old['d_full']-(old['d_steer']+old['d_loc']+old['d_life'])

    rows[0]['summed_single_factor_residual_tv'] = float(abs(resid).sum()/2)

    rows[1]['summed_single_factor_residual_tv'] = float(abs(resid_old).sum()/2)

    pd.DataFrame(rows).to_csv(_LAYOUT_ROOT / 'results/models/bam_primary_calibration_fit.csv', index=False)

    print(pd.DataFrame(rows).to_string())



if __name__ == '__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()

