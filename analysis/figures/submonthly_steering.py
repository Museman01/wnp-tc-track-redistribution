"""Draw submonthly steering regimes and track-composition diagnostics.

Usage: python -m analysis.figures.submonthly_steering <project_root> <gshhg_i_L1.shp> <out_dir>
Style follows the SI map figures: one font size for ordinary text, bold panel letters,
GSHHG intermediate L1 land, RdBu_r diverging colours, no drawn figure title.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
import geopandas as gpd

ROOT, SHP, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
names = {f.name for f in font_manager.fontManager.ttflist}
FONT = next((f for f in ['Arial', 'Liberation Sans', 'DejaVu Sans'] if f in names), 'DejaVu Sans')
SIZE = 13
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': [FONT, 'DejaVu Sans'], 'font.size': SIZE,
                     'axes.titlesize': SIZE, 'axes.labelsize': SIZE, 'xtick.labelsize': SIZE, 'ytick.labelsize': SIZE,
                     'legend.fontsize': SIZE, 'mathtext.fontset': 'custom', 'mathtext.rm': FONT, 'mathtext.it': FONT + ':italic',
                     'axes.linewidth': .8, 'pdf.fonttype': 42, 'ps.fonttype': 42})
LAND = gpd.read_file(SHP, bbox=(80, -5, 185, 60))
REGIME_COLORS = ['#355C7D', '#6C9A3B', '#B5523B', '#8E6BB0']

def geo(ax, extent):
    LAND.plot(ax=ax, facecolor='#E6E4DF', edgecolor='#505050', linewidth=.35, zorder=3)
    ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel(''); ax.set_ylabel('')
    ax.set_aspect(1 / np.cos(np.deg2rad((extent[2] + extent[3]) / 2)), adjustable='box')
    ax.set_xticks(range(100, 181, 20)); ax.set_yticks(range(0, extent[3] + 1, 10))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: '180°' if np.isclose(v, 180) else f'{v:g}°E'))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: '0°' if np.isclose(v, 0) else f'{v:g}°N'))

def letter(ax, s):
    ax.text(.015, .98, s, transform=ax.transAxes, ha='left', va='top', fontsize=SIZE, fontweight='bold', zorder=12,
            bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=.6))

def save(fig, name):
    fig.savefig(OUT / f'{name}.png', dpi=300)
    fig.savefig(OUT / f'{name}.pdf', metadata={'Creator': 'Matplotlib', 'Author': 'Daizhao Huang; Yuanfu Yue'})
    plt.close(fig)

# ------------------------------------------------------------------ Figure S20
a = np.load(ROOT / 'data/figure_inputs/submonthly_steering/figS20_regimes.npz')
lon, lat, comp = a['longitude'], a['latitude'], a['composite_uv_anomaly']
fig = plt.figure(figsize=(3332.01 / 300, 2860.01 / 300))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, .80], left=.075, right=.86, bottom=.07, top=.985, hspace=.24, wspace=.16)
lim = 8.0
for k in range(4):
    ax = fig.add_subplot(gs[k // 2, k % 2])
    im = ax.pcolormesh(np.r_[lon - 1.25, lon[-1] + 1.25], np.r_[lat - 1.25, lat[-1] + 1.25], comp[k, 0], cmap='RdBu_r',
                       vmin=-lim, vmax=lim, shading='flat', zorder=1)
    s = slice(None, None, 2)
    q = ax.quiver(lon[s], lat[s], comp[k, 0][s, s], comp[k, 1][s, s], scale=60, width=.0035, color='#1A1A1A', zorder=5)
    geo(ax, (100, 180, 0, 45)); letter(ax, 'abcd'[k])
    ax.text(.985, .035, f'Regime {k + 1}\n1966–1995: {a["early_frequency"][k]:.3f}\n1996–2025: {a["late_frequency"][k]:.3f}',
            transform=ax.transAxes, ha='right', va='bottom', zorder=12, linespacing=1.15,
            bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=2.5))
    if k == 0:
        qk = ax.quiverkey(q, .30, .065, 4, r'4 m s$^{-1}$', labelpos='E', coordinates='axes', zorder=13)
        qk.text.set_bbox(dict(facecolor='white', edgecolor='none', alpha=.9, pad=1.5))
cax = fig.add_axes([.885, .42, .018, .52])
fig.colorbar(im, cax=cax, label='Zonal-wind anomaly (m s$^{-1}$)', extend='both')
ax = fig.add_subplot(gs[2, :])
years = a['years']
for k in range(4):
    ax.plot(years, a['annual_frequency'][:, k], color=REGIME_COLORS[k], lw=1.4, label=f'Regime {k + 1}')
ax.axvline(1995.5, color='#666666', lw=.9, ls=':')
ax.set_xlim(1965.5, 2025.5); ax.set_ylabel('Regime frequency'); ax.set_xlabel('Year')
ax.spines[['top', 'right']].set_visible(False)
ax.legend(ncol=4, loc='upper left', frameon=False, bbox_to_anchor=(0, 1.13), handlelength=1.6, columnspacing=1.2)
letter(ax, 'e')
save(fig, 'FigS20_steering_regimes')

# ------------------------------------------------------------------ Figure S21
b = np.load(ROOT / 'data/figure_inputs/submonthly_steering/figS21_decomposition.npz')
ks = pd.read_csv(ROOT / 'results/circulation/submonthly_steering/k_sensitivity_path_decomposition.csv')
proj = {'frequency': (0.049483, -2.936767, 3.606184), 'within': (99.950517, 81.174067, 118.690094)}
fig = plt.figure(figsize=(3332.01 / 300, 1980.01 / 300))
gs = fig.add_gridspec(2, 2, left=.06, right=.985, bottom=.235, top=.975, hspace=.28, wspace=.22)
lim = 0.70
maps = {}
for i, (key, pos) in enumerate([('observed_pp', (0, 0)), ('frequency_pp', (0, 1)), ('within_pp', (1, 0))]):
    ax = fig.add_subplot(gs[pos])
    cmap, l = ('PuOr_r', lim / 20) if key == 'frequency_pp' else ('RdBu_r', lim)
    maps[key] = (ax, ax.pcolormesh(b['lon_edges'], b['lat_edges'], b[key], cmap=cmap, vmin=-l, vmax=l, zorder=1))
    geo(ax, (100, 180, 0, 40)); letter(ax, 'abc'[i])
ax = fig.add_subplot(gs[1, 1])
x = np.array([0, 1]); vals = [proj['frequency'], proj['within']]
for xi, (v, lo, hi), c in zip(x, vals, ['#6C7A89', '#2F4B7C']):
    ax.bar(xi, v, width=.55, color=c, alpha=.85, zorder=2)
    ax.errorbar(xi, v, yerr=[[v - lo], [hi - v]], fmt='none', ecolor='#1A1A1A', elinewidth=1.1, capsize=5, zorder=3)
refit = ks[ks.labels.str.startswith('refit')]
off = np.linspace(-.18, .18, len(refit))
ax.scatter(0 + off, refit.frequency_percent, s=34, facecolor='white', edgecolor='#B5523B', lw=1.3, zorder=4)
ax.axhline(0, color='#555555', lw=.8)
ax.set_xticks(x); ax.set_xticklabels(['Regime\nfrequency', 'Within\nregime'])
ax.set_ylabel('Share of projected contrast (%)'); ax.set_ylim(-15, 130); ax.set_xlim(-.6, 1.6)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(handles=[Line2D([0], [0], marker='o', ls='none', mfc='white', mec='#B5523B', mew=1.3, ms=6, label='Refit, k = 4–8')],
          loc='upper left', frameon=False, bbox_to_anchor=(0.10, 1.0))
letter(ax, 'd')
ax_d = ax
fig.canvas.draw()
for host, key, label, ticks in [(maps['within_pp'][0], 'within_pp', 'Panels a and c (percentage points)', [-.6, -.3, 0, .3, .6]),
                                (ax_d, 'frequency_pp', 'Panel b (percentage points)', [-.03, -.02, -.01, 0, .01, .02, .03])]:
    bb = host.get_position()
    cax = fig.add_axes([bb.x0, .105, bb.width, .024])
    cb = fig.colorbar(maps[key][1], cax=cax, orientation='horizontal', extend='both', label=label)
    cb.set_ticks(ticks); cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _, n=(1 if abs(ticks[-1]) > .1 else 2): f'{v:.{n}f}'.replace('-', '\u2212')))
save(fig, 'FigS21_regime_track_decomposition')
print('font', FONT)
