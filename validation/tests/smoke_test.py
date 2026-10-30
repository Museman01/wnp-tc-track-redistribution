"""用极小合成数据跑通活动 Core 主链路并检查基本不变量。

运行：python -m validation.tests.smoke_test  （或 pytest validation/tests/smoke_test.py）
"""
import numpy as np
import pandas as pd


def test_smoke_core_chain():
    """Core 主链路冒烟：玩具轨迹上的登陆几何与统计工具应正常工作、量级合理。"""
    from core.landfall_gshhg import densify, is_land
    from core.stats_utils import effective_n, theil_sen_ci
    from core.wavelet import cwt_morlet, global_spectrum

    # 1) 登陆几何：海->陆直线恰一次穿越。
    res = 0.02
    lon0, lat1 = 120.0, 20.0
    nx, ny = int((140 - 120) / res), int((20 - 10) / res)
    mask = np.zeros((ny, nx), dtype="uint8")
    mask[:, int((130 - lon0) / res):] = 1
    track = pd.DataFrame({"lon": [125., 131.], "lat": [15., 15.]})
    f, lat, lon = densify(track.iloc[0], track.iloc[1], res / 2)
    state = [is_land(mask, lon0, lat1, res, la, lo) for la, lo in zip(lat, lon)]
    crossings = sum(1 for j in range(1, len(state)) if state[j] and not state[j - 1])
    assert crossings == 1

    # 2) 统计工具在合成年序列上可运行且产出合理。
    yrs = np.arange(1980, 2025)
    series = 0.3 * (yrs - 1980) + np.sin(np.arange(len(yrs)))
    assert effective_n(series, series) >= 3
    slope, lo, hi = theil_sen_ci(yrs, series, nboot=200, seed=0)
    assert lo <= slope <= hi

    # 3) 小波在合成 4 周期信号上识别出 ~4 的峰值周期。
    sig = np.sin(2 * np.pi * np.arange(256) / 4.0)
    power, period, scale, coi = cwt_morlet(sig, dt=1.0)
    assert abs(period[int(np.argmax(global_spectrum(power)))] - 4.0) < 1.0


if __name__ == "__main__":
    test_smoke_core_chain()
    print("core smoke test passed")
