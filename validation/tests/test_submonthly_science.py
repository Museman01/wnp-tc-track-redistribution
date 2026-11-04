import numpy as np
import pandas as pd
from analysis.circulation import submonthly_statistics as s
from analysis.circulation.sampling_sensitivity import format_s7_pvalue


def test_symmetric_path_partition_closes():
    early=np.array([[[0.1,0.2]],[[0.3,0.4]]])
    late=np.array([[[0.3,0.1]],[[0.2,0.4]]])
    result=s.decompose(early,late)
    np.testing.assert_allclose(result['between']+result['within'],result['delta'],atol=1e-15)
    assert abs(s.project(result['between'],result['delta'])+s.project(result['within'],result['delta'])-100)<1e-12


def test_incomplete_environmental_annulus_has_no_partial_mean():
    latitude=np.arange(60,-11,-1);longitude=np.arange(90,191)
    weights,coverage=s.geometry_weights(latitude,longitude,15,180,'annulus_5_7')
    assert coverage['complete_geometry'] and np.isclose(weights.sum(),1)
    weights,coverage=s.geometry_weights(latitude,longitude,15,189,'annulus_5_7')
    assert weights is None and not coverage['complete_geometry']


def test_motion_segment_uses_short_longitude_difference():
    frame=pd.DataFrame(dict(sid=['one','one'],season=[2000,2000],iso_time=pd.to_datetime(['2000-06-01','2000-06-01 06:00'],format='mixed'),lat=[10.,10.],lon=[179.,-179.]))
    result=s.motion_pairs(frame)
    assert len(result)==1 and result.obs_u.iloc[0]>0 and result.obs_v.iloc[0]==0


def test_s7_decimal_display_preserves_half_up_boundary():
    assert format_s7_pvalue('0.0825')=='0.083'
    assert format_s7_pvalue('0.0075')=='0.008'
    assert format_s7_pvalue('0.1795')=='0.180'
    assert format_s7_pvalue('0.0001')=='<0.001'
