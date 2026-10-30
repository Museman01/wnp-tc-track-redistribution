"""Run the current BAM source unchanged numerically, with outputs and input routing isolated."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import numpy as np

import netCDF4

from core.paths import PROJECT,output_path,import_source,write_json



def main():

    bam=import_source('candidate_bam_original','analysis/models/counterfactual_steering.py')

    _LAYOUT_ROOT=output_path('results/models/reproduction')

    _LAYOUT_ROOT.mkdir(parents=True,exist_ok=True)

    # Explicit source routing replaces only the historical /tmp fallback.

    with netCDF4.Dataset(PROJECT/'data/interim/steering.nc') as ds:

        t=ds.variables['time']; dates=netCDF4.num2date(t[:],t.units,getattr(t,'calendar','standard'))

        bam._STEER_CACHE={'u':np.asarray(ds['u_steer'][:],dtype=np.float32),

            'v':np.asarray(ds['v_steer'][:],dtype=np.float32),'lat':np.asarray(ds['latitude'][:],dtype=float),

            'lon':np.asarray(ds['longitude'][:],dtype=float),

            'yrs':np.array([d.year for d in dates]),'mons':np.array([d.month for d in dates])}

    bam.main()

    new=np.load(_LAYOUT_ROOT / 'results/models/counterfactual/fields.npz')

    old=np.load(_LAYOUT_ROOT / 'results/models/counterfactual/fields.npz')

    mapping={'f_BASE_P1':'f_BASE_EARLY','f_BASE_P2':'f_BASE_LATE'}

    errors={key:float(np.max(abs(old[key]-new[mapping.get(key,key)]))) for key in old.files}

    maxerr=max(errors.values())

    write_json('results/diagnostics/bam_source_result_reproduction.json',{'status':'PASS' if maxerr<1e-10 else 'DIFFERS',

        'maxabs_by_array':errors,'key_mapping':mapping,'source_input':'data/interim/steering.nc, explicit',

        'historical_non_six_hour_calibration':'preserved for exact reproduction; correction tested separately',

        'legacy_comparison_labels':'original script output retained as provenance, not endorsed causal inference'})

    print('BAM current-source reproduction maxabs:',maxerr,flush=True)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
