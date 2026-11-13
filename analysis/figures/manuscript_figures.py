"""Render the five data figures used in the manuscript."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Soft teal-white-coral diverging map for the within-regime term (Figure 4f); zero is pure white.
TEAL_CORAL = LinearSegmentedColormap.from_list('teal_coral', ['#1B6F73', '#7FBFBE', '#FFFFFF', '#F2A68C', '#B5452F'])
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
from core.paths import PROJECT, output_path
from analysis.figures.main_figures import (COLORS, DISPLAY, AG3, AG4, ZA, ZV, H, add_coast, format_geo_axes,
                                           open_panel_axis, scale_typography)
from analysis.figures.shared_figure_style import PERIODS

OUT = None
MODE = 'cells'   # Figure 5 maps: 1-degree share changes ('points' draws individual crossings)
EARLY, LATE = PERIODS['early'][0], PERIODS['late'][0]


def letter(ax, s):
    ax.text(0.015, 0.985, s, transform=ax.transAxes, va='top', ha='left', fontsize=12, fontweight='bold', zorder=12)


def geo(ax, extent=(100, 0, 180, 40), last=True):
    add_coast(ax, extent); format_geo_axes(ax)
    ax.set_xticks(list(range(int(extent[0]), int(extent[2]) + (1 if last else 0), 20)))
    ax.set_yticks(list(range(int(extent[1]), int(extent[3]) + 1, 10 if extent[3] - extent[1] <= 40 else 20)))


def save(fig, stem, number, title, scale=1.0, dpi=300):
    for ext in ('.png', '.pdf'):
        if (OUT / (stem + ext)).exists(): raise FileExistsError(stem)
    scale_typography(fig, scale)
    fig.savefig(OUT / f'{stem}.png', dpi=dpi, bbox_inches='tight')
    fig.savefig(OUT / f'{stem}.pdf', bbox_inches='tight',
                metadata={'Title': title, 'Subject': f'Journal of Geophysical Research: Atmospheres Figure {number}', 'Creator': 'Matplotlib'})
    plt.close(fig); print('Figure', stem, flush=True)


# ---------------------------------------------------------------- Figure 1: series (a, b) and four agency maps (c-f)
def fig1():
    annual = pd.read_csv(ZA / '03_common_storms/core_crossagency_annual.csv'); hf = H('frequency.csv')
    idx = pd.read_csv(ZV / 'supplemental_validation/results/wnp_tc_redistribution_index_annual.csv'); hi = H('index.csv')
    tr = np.load(PROJECT / 'results/tracks/hko_track_fields.npz')
    fig = plt.figure(figsize=(11.6, 7.9))
    gs = fig.add_gridspec(3, 3, width_ratios=[1, 1, 0.035], height_ratios=[0.95, 1, 1], left=0.07, right=0.93,
                          bottom=0.05, top=0.98, hspace=0.30, wspace=0.10)
    axa = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1])
    for a in AG3:
        p = annual[annual.agency.eq(a)]; axa.plot(p.season, p.n_tc, color=COLORS[a], lw=1.2, label=DISPLAY[a])
    axa.plot(hf.season, hf.n_tc, color=COLORS['HKO'], lw=1.2, label='HKO')
    axa.set(xlabel='Year', ylabel='TC count'); axa.legend(frameon=False, ncol=4, loc='upper right')
    axa.set_ylim(top=axa.get_ylim()[1] + 6)
    for a in ['PRIMARY'] + AG3:
        d = idx[(idx.agency == a) & (idx.weighting == 'track_point')].sort_values('year')
        axb.plot(d.year, d.index_oos, lw=1.5 if a == 'PRIMARY' else 1.1, color=COLORS[a], label=DISPLAY[a], zorder=4 if a == 'PRIMARY' else 3)
    axb.plot(hi.year, hi.index_oos, lw=1.1, color=COLORS['HKO'], label='HKO')
    axb.axhline(0, color='0.45', lw=0.7); axb.axvline(1995.5, color='0.55', ls='--', lw=0.8)
    axb.set(xlabel='Year', ylabel='Leave-one-year-out\nredistribution index'); axb.legend(frameon=False, ncol=3, loc='upper left', bbox_to_anchor=(0.05, 1.0))
    axb.set_ylim(top=axb.get_ylim()[1] + 0.25 * np.ptp(axb.get_ylim()))
    open_panel_axis(axa); open_panel_axis(axb, right=True)
    lim = max(abs(tr['native_all_' + k]).max() for k in AG4) * 100
    maps = []
    for i, k in enumerate(AG4):
        ax = fig.add_subplot(gs[1 + i // 2, i % 2])
        m = ax.imshow(tr['native_all_' + k] * 100, origin='lower', extent=[100, 180, 0, 40], cmap='RdBu_r', vmin=-lim, vmax=lim, aspect='auto', zorder=1)
        geo(ax, last=bool(i % 2)); maps.append(ax)
        if i % 2: ax.tick_params(axis='y', labelleft=False)
        if i < 2: ax.tick_params(axis='x', labelbottom=False)
    for ax, s in zip([axa, axb] + maps, 'abcdef'): letter(ax, s)
    fig.canvas.draw()
    top = maps[1].get_position(); bot = maps[3].get_position()
    cax = fig.add_axes([top.x1 + 0.012, bot.y0, 0.014, top.y1 - bot.y0])
    fig.colorbar(m, cax=cax, label='Change (percentage points)')
    save(fig, 'Fig01_cross_agency_redistribution', 1, 'Cross-agency frequency, redistribution index and track-composition change', 1.70)


# ---------------------------------------------------------------- Figure 2: genesis-distribution and conditional-track-field terms
def fig2():
    f = np.load(PROJECT / 'results/genesis/PRIMARY_fields.npz'); lo, la = f['lon_edges'], f['lat_edges']
    arr = {k: f[f'{k}_field'].reshape(len(la) - 1, len(lo) - 1) * 100 for k in ('genesis', 'propagation')}
    fig = plt.figure(figsize=(11.6, 3.75))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.06], left=0.06, right=0.985, bottom=0.10, top=0.98, hspace=0.40, wspace=0.07)
    axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]; ms = []
    for ax, k, name, s, cmap in zip(axs, ['genesis', 'propagation'], ['Genesis distribution', 'Conditional track field'], 'ab', ['PuOr_r', 'RdBu_r']):
        v = np.abs(arr[k]).max()
        ms.append(ax.pcolormesh(lo, la, arr[k], cmap=cmap, vmin=-v, vmax=v, shading='flat', zorder=1))
        geo(ax, last=(s == 'b')); letter(ax, s)
    axs[1].tick_params(axis='y', labelleft=False)
    fig.canvas.draw()
    for ax, mm, lab in [(axs[0], ms[0], 'Contribution (percentage points per cell)'), (axs[1], ms[1], 'Contribution (percentage points per cell)')]:
        p = ax.get_position(); cax = fig.add_axes([p.x0 + 0.1 * p.width, 0.075, 0.8 * p.width, 0.035])
        cb = fig.colorbar(mm, cax=cax, orientation='horizontal', label=lab)
        cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.2f}'.rstrip('0').rstrip('.').replace('-', '\u2212') if v else '0'))
    save(fig, 'Fig02_genesis_conditional_track_field_terms', 2, 'Genesis-distribution and conditional-track-field terms', 1.70, 320)


# ---------------------------------------------------------------- Figure 3: index regressions (a, b) and two-period circulation (c, d)
def fig3():
    F = np.load(ZV / 'circulation_metric_validation/results/wnp_tc_circulation_fields.npz')
    S = np.load(PROJECT / 'results/circulation/steering_field_regressions.npz')
    lon, lat = F['longitude'], F['latitude']; X, Y = np.meshgrid(S['longitude'], S['latitude'])
    vmax = max(np.abs(F['raw_z_eddy_beta']).max(), np.abs(F['detrended_z_eddy_beta']).max())
    zlim = float(np.nanmax(abs(S['period_change_z_eddy'])))
    fig = plt.figure(figsize=(11.6, 7.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.035], left=0.07, right=0.92, bottom=0.05, top=0.98, hspace=0.06, wspace=0.08)
    ax = [[fig.add_subplot(gs[r, c]) for c in range(2)] for r in range(2)]
    for a_, sc, name, s in [(ax[0][0], 'raw', 'Raw regression', 'a'), (ax[0][1], 'detrended', 'Detrended regression', 'b')]:
        m1 = a_.pcolormesh(lon, lat, F[f'{sc}_z_eddy_beta'], cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto', zorder=1)
        # Reference length (stated in the caption): 1 m s-1 per index SD spans 6 degrees of longitude.
        a_.quiver(lon[::2], lat[::2], F[f'{sc}_u_beta'][::2, ::2], F[f'{sc}_v_beta'][::2, ::2], color='0.15',
                  scale=1 / 6, scale_units='x', width=0.0024, zorder=5)
        geo(a_, (80, 0, 180, 65), last=(s == 'b')); letter(a_, s)
    # Panel c (zero-wind lines and steering change) on the left; panel d (height change) beside its colour bar.
    d, c = ax[1]
    m2 = c.pcolormesh(S['longitude'], S['latitude'], S['period_change_z_eddy'], cmap='PuOr_r', vmin=-zlim, vmax=zlim, shading='auto', zorder=1)
    for p, col in [('early', EARLY), ('late', LATE)]:
        cs = c.contour(S['longitude'], S['latitude'], S['z500_' + p], levels=[5880], colors=[col], linewidths=1.3, zorder=6)
        d.contour(S['longitude'], S['latitude'], S['u500_' + p], levels=[0], colors=[col], linewidths=1.3, zorder=6)
    # Reference length (stated in the caption, panel c): 1 m s-1 spans 7.5 degrees of longitude.
    d.quiver(X[::2, ::2], Y[::2, ::2], S['period_change_u300'][::2, ::2], S['period_change_v300'][::2, ::2],
             scale=1 / 7.5, scale_units='x', width=0.0027, color='0.15', zorder=5)
    geo(d, (80, 0, 180, 65), last=False); geo(c, (80, 0, 180, 65))
    letter(d, 'c'); letter(c, 'd')
    for a_ in (ax[0][1], c): a_.tick_params(axis='y', labelleft=False)
    for a_ in ax[0]: a_.tick_params(axis='x', labelbottom=False)
    fig.canvas.draw()
    for row, mm, lab in [(0, m1, 'Regression (m per index SD)'),
                         (1, m2, 'Change (m)')]:
        pb = ax[row][1].get_position(); cax = fig.add_axes([pb.x1 + 0.012, pb.y0, 0.014, pb.height])
        fig.colorbar(mm, cax=cax, label=lab)
    save(fig, 'Fig03_circulation_regression_and_two_period_change', 3, 'Circulation regressions and two-period circulation change', 1.68)


# ---------------------------------------------------------------- Figure 4: regime-frequency and within-regime terms
def fig4():
    """Regime composites (a-d) and the regime decomposition of the June-October track-composition difference (e, f)."""
    A = np.load(PROJECT / 'data/figure_inputs/submonthly_steering/figS20_regimes.npz')
    b = np.load(PROJECT / 'data/figure_inputs/submonthly_steering/figS21_decomposition.npz'); lim = 0.75   # covers the largest within-regime value (0.716) without colour-bar extensions
    lon, lat, comp = A['longitude'], A['latitude'], A['composite_uv_anomaly']
    fig = plt.figure(figsize=(11.6, 10.4))
    gs = fig.add_gridspec(4, 3, width_ratios=[1, 1, 0.035], height_ratios=[1, 1, 0.93, 0.05], left=0.06, right=0.93,
                          bottom=0.05, top=0.985, hspace=0.16, wspace=0.07)
    for k in range(4):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        im = ax.pcolormesh(np.r_[lon - 1.25, lon[-1] + 1.25], np.r_[lat - 1.25, lat[-1] + 1.25], comp[k, 0], cmap='RdBu_r',
                           vmin=-9, vmax=9, shading='flat', zorder=1)
        sl = slice(None, None, 2)
        # Reference length (stated in the caption): 1 m s-1 spans 1.5 degrees of longitude.
        ax.quiver(lon[sl], lat[sl], comp[k, 0][sl, sl], comp[k, 1][sl, sl], scale=1 / 1.5, scale_units='x', width=0.0032, color='#1A1A1A', zorder=5)
        geo(ax, (100, 0, 180, 45), last=bool(k % 2)); letter(ax, 'abcd'[k])
        if k % 2: ax.tick_params(axis='y', labelleft=False)
        ax.tick_params(axis='x', labelbottom=False)
    fig.canvas.draw()
    top = fig.axes[1].get_position(); bot = fig.axes[3].get_position()
    cb = fig.colorbar(im, cax=fig.add_axes([top.x1 + 0.012, bot.y0, 0.014, top.y1 - bot.y0]))
    cb.set_label('Zonal-wind anomaly (m s$^{-1}$)')
    axs = [fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])]; ms = []
    for ax, key, name, s in zip(axs, ['frequency_pp', 'within_pp'], ['Regime-frequency term', 'Within-regime term'], 'ef'):
        cmap, l = ('PuOr_r', lim / 20) if key == 'frequency_pp' else (TEAL_CORAL, lim)
        ms.append(ax.pcolormesh(b['lon_edges'], b['lat_edges'], b[key], cmap=cmap, vmin=-l, vmax=l, zorder=1))
        geo(ax, last=(s == 'f')); letter(ax, s)
    axs[1].tick_params(axis='y', labelleft=False)
    fig.canvas.draw()
    for ax, mm, ticks, lab in [(axs[0], ms[0], [-.03, 0, .03], 'Contribution (percentage points)'),
                               (axs[1], ms[1], [-.6, -.3, 0, .3, .6], 'Contribution (percentage points)')]:
        p = ax.get_position(); cax = fig.add_axes([p.x0 + 0.1 * p.width, p.y0 - 0.068, 0.8 * p.width, 0.014])
        cb = fig.colorbar(mm, cax=cax, orientation='horizontal', label=lab); cb.set_ticks(ticks)
        n = 2 if abs(ticks[-1]) < .1 else 1
        cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _, n=n: f'{v:.{n}f}'.replace('-', '\u2212')))
    save(fig, 'Fig04_steering_regimes_and_track_decomposition', 4, 'Steering regimes and regime decomposition of track composition', 1.70, 320)


# ---------------------------------------------------------------- Figure 5: landfall-location maps and latitude distributions, (a) all events and (b) first landfall
TOPO = PROJECT / 'data/derived/etopo1_99-146E_4-46N_2min.nc'   # ETOPO1 subset (Amante & Eakins, 2009)


def landfall_events():
    e = pd.read_csv(PROJECT / 'results/landfall/landfall_events_exact.csv')
    h = pd.read_csv(PROJECT / 'results/landfall/hko_landfall_events.csv')
    cols = ['agency', 'sid', 'season', 'time', 'lat', 'lon', 'pre_crossing_native_ts']
    ev = pd.concat([e[cols], h[cols]])
    ev = ev[ev.pre_crossing_native_ts.astype(str).eq('True') & ev.season.between(1966, 2025)].copy()
    ev['period'] = np.where(ev.season <= 1995, 'early', 'late')
    first = ev.sort_values('time').groupby(['agency', 'sid']).head(1)
    ann = pd.read_csv(output_path('results/agency_comparison/landfall/landfall_annual_four_agency.csv'))
    for d, x in [('all_events', ev), ('first_landfall', first)]:   # events reproduce the annual means in the text
        m = x.groupby(['agency', 'season']).lat.mean().reset_index()
        a = ann[ann.definition.eq(d)].merge(m, left_on=['agency', 'year'], right_on=['agency', 'season'])
        assert len(a) == 240 and float((a.mean_lat - a.lat).abs().max()) < 1e-9, d
    return ev, first, ann


def relief(ax, extent):
    import xarray as xr, geopandas as gpd
    from matplotlib.colors import LightSource, LinearSegmentedColormap
    t = xr.open_dataset(TOPO).altitude.sel(latitude=slice(extent[1], extent[3]), longitude=slice(extent[0], extent[2]))
    z = t.values.astype(float); z = z[::-1] if t.latitude.values[0] < t.latitude.values[-1] else z
    ls = LightSource(azdeg=315, altdeg=40)
    land = LinearSegmentedColormap.from_list('land', ['#C9D3AE', '#D9D2A6', '#C7B48C', '#A78F72', '#8C7A68', '#F1EEEA'])
    sea = LinearSegmentedColormap.from_list('sea', ['#9FB7CB', '#BFD0DE', '#DCE7EF', '#EEF3F7'])
    rgb_land = ls.shade(np.clip(z, 0, 4500), cmap=land, vmin=0, vmax=4500, vert_exag=0.04, blend_mode='soft')
    rgb_sea = ls.shade(np.clip(z, -7000, 0), cmap=sea, vmin=-7000, vmax=0, vert_exag=0.006, blend_mode='soft')
    # Muted relief: blend the shaded colours towards white so tracks and landfall points remain dominant.
    rgb = np.where((z > 0)[..., None], 0.45 * rgb_land + 0.55, 0.6 * rgb_sea + 0.4)
    ax.imshow(rgb, extent=[float(t.longitude.min()), float(t.longitude.max()), float(t.latitude.min()), float(t.latitude.max())],
              origin='upper', interpolation='bilinear', zorder=1)
    coast = gpd.read_file(PROJECT / 'data/raw/GSHHG/GSHHS_shp/i/GSHHS_i_L1.shp', bbox=(extent[0], extent[1], extent[2], extent[3]))
    coast.boundary.plot(ax=ax, color='#4A4A4A', linewidth=0.35, zorder=3)
    ax.set_xlim(extent[0], extent[2]); ax.set_ylim(extent[1], extent[3])
    ax.set_aspect(1 / np.cos(np.deg2rad((extent[1] + extent[3]) / 2)), adjustable='box'); format_geo_axes(ax)


def fig5():
    from scipy.stats import gaussian_kde
    ev, first, ann = landfall_events(); ext = (100, 5, 145, 45); y = np.linspace(ext[1], ext[3], 801)
    fig = plt.figure(figsize=(9.6, 10.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 0.62], left=0.08, right=0.97, bottom=0.075, top=0.985, hspace=0.07, wspace=0.035)
    rng = np.random.default_rng(1996)
    from matplotlib.collections import LineCollection
    tr = pd.read_csv(PROJECT / 'data/raw/ibtracs_wp.csv', usecols=['SID', 'SEASON', 'ISO_TIME', 'LAT', 'LON', 'TRACK_TYPE'], skiprows=[1], low_memory=False)
    tr = tr[tr.TRACK_TYPE.eq('main') & tr.SID.isin(set(ev.sid))].copy()
    tr['LAT'] = pd.to_numeric(tr.LAT); tr['LON'] = pd.to_numeric(tr.LON); tr['SEASON'] = pd.to_numeric(tr.SEASON)
    segs, cols = [], []
    for sid, g in tr.sort_values(['SID', 'ISO_TIME']).groupby('SID'):   # IBTrACS combined positions, one track per landfalling storm
        segs.append(g[['LON', 'LAT']].to_numpy()); cols.append(EARLY if g.SEASON.iloc[0] <= 1995 else LATE)
    order = rng.permutation(len(segs)); segs = [segs[k] for k in order]; cols = [cols[k] for k in order]
    print('landfalling storms with tracks:', len(segs))
    for r, (x, d, name) in enumerate([(ev, 'all_events', 'All landfall events'), (first, 'first_landfall', 'First landfall')]):
        axm = fig.add_subplot(gs[r, 0]); axp = fig.add_subplot(gs[r, 1], sharey=axm)
        relief(axm, ext)
        axm.add_collection(LineCollection(segs, colors=cols, linewidths=0.5, alpha=0.10, zorder=4))
        if MODE == 'points':
            q = x.iloc[rng.permutation(len(x))]   # random drawing order so neither period hides the other
            axm.scatter(q.lon, q.lat, s=11, c=np.where(q.period.eq('early'), EARLY, LATE), alpha=0.8, edgecolors='white', linewidths=0.25, zorder=5)
        else:
            # Share of each period's landfalls in 1-degree cells, averaged over the four records; circle area ~ |late - early|.
            x = x.assign(ci=np.floor(x.lat).astype(int), cj=np.floor(x.lon).astype(int))
            sh = []
            for a in AG4:
                g = x[x.agency.eq(a)]
                t = g.groupby(['period', 'ci', 'cj']).size().unstack('period', fill_value=0)
                t = t / t.sum(); sh.append(t)
            t = pd.concat(sh).groupby(level=[0, 1]).sum() / len(AG4); t = t.reindex(columns=['early', 'late'], fill_value=0)
            dlt = (t.late - t.early) * 100; dlt = dlt[dlt.abs() > 0.05]
            yy = dlt.index.get_level_values(0) + 0.5; xx = dlt.index.get_level_values(1) + 0.5
            axm.scatter(xx, yy, s=np.abs(dlt.values) * 160, c=np.where(dlt.values > 0, LATE, EARLY), alpha=0.75,
                        edgecolors='white', linewidths=0.5, zorder=5)
            if r == 0: print('max |change| (pp per cell):', float(dlt.abs().max()))
        dens = {}
        for p in ('early', 'late'):
            sub = x[x.period.eq(p)]
            dens[p] = np.mean([gaussian_kde(sub[sub.agency.eq(a)].lat, bw_method=0.18)(y) for a in AG4], axis=0) * 100
        e_, l_ = dens['early'], dens['late']
        axp.fill_betweenx(y, e_, l_, where=l_ >= e_, color=LATE, alpha=0.45, lw=0, interpolate=True, zorder=2)
        axp.fill_betweenx(y, e_, l_, where=l_ < e_, color=EARLY, alpha=0.45, lw=0, interpolate=True, zorder=2)
        axp.fill_betweenx(y, 0, np.minimum(e_, l_), color='0.88', lw=0, zorder=1)
        axp.plot(e_, y, color=EARLY, lw=1.6, zorder=3); axp.plot(l_, y, color=LATE, lw=1.6, zorder=3)
        for p, c, yr in [('early', EARLY, (1966, 1995)), ('late', LATE, (1996, 2025))]:
            mu = np.mean([ann[ann.definition.eq(d) & ann.agency.eq(a) & ann.year.between(*yr)].mean_lat.mean() for a in AG4])
            axp.axhline(mu, color=c, ls='--', lw=1.2, zorder=4)
        for v in (10, 20, 30, 40): axp.axhline(v, color='0.85', lw=0.6, zorder=0)
        axp.set_xlim(0, max(e_.max(), l_.max()) * 1.08); axp.tick_params(axis='y', labelleft=False)
        axp.spines[['top', 'right']].set_visible(False)
        axm.set_xticks([100, 110, 120, 130, 140]); axm.set_yticks([10, 20, 30, 40])
        if r == 0: axm.tick_params(axis='x', labelbottom=False); axp.tick_params(axis='x', labelbottom=False)
        else: axp.set_xlabel('Share of landfalls (% per degree latitude)')
        letter(axm, 'ab'[r])   # map and latitude profile of one landfall definition form one panel
    save(fig, 'Fig05_landfall_share_change_and_latitude_distribution', 5, 'Landfall locations and latitude distributions', 1.40)


def main():
    import argparse
    global OUT, MODE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('figures', nargs='*', choices=['1', '2', '3', '4', '5'])
    parser.add_argument('--points', action='store_true')
    args = parser.parse_args()
    if not args.execute: parser.error('add --execute to render figures')
    OUT = args.output_dir.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.points: MODE = 'points'
    for number in (args.figures or ['1', '2', '3', '4', '5']): globals()['fig'+number]()


if __name__ == '__main__':
    main()
