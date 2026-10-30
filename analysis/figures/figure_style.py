"""AGU 版式与全工程统一绘图风格。

所有活动出图函数复用本模块，保持配色、字体和版面一致。
"""

from core.paths import PROJECT as _LAYOUT_ROOT
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from core.utils import load_config, ROOT

# ---- 字体与版式（AGU 图件统一使用 Arial 无衬线字体）----
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Microsoft YaHei"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial",
    "axes.unicode_minus": False,
    # The final figures are commonly read at manuscript-page scale.  The
    # previous 7--9 pt settings were suitable for screen previews but made
    # panel labels and axis text undersized after placement in Word/PDF.
    "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
    "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "legend.fontsize": 9.5,
    "axes.titlepad": 7,
    "figure.titlesize": 13, "figure.titleweight": "bold",
    "axes.linewidth": 0.75, "xtick.major.width": 0.75, "ytick.major.width": 0.75,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "figure.dpi": 150, "savefig.dpi": 600,
    "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,  # 可编辑矢量字体
})

# ---- 期刊版面宽度（英寸）----
COL1, COL1_5, COL2 = 3.5, 5.0, 7.2     # 单栏 89mm / 1.5 栏 / 双栏 183mm

# ---- 莫奈分类色（离散用：分型、位相、折线、柱状）----
# Use a conventional academic blue for the primary series; the remaining
# muted categorical colours retain the existing manuscript palette.
MONET = ["#4C72B0", "#7E9CC4", "#9E8FC0", "#8FB996",
         "#E0908E", "#E0B97D", "#B57BA6", "#8A93A6"]
mpl.rcParams["axes.prop_cycle"] = mpl.cycler(color=MONET)

# ---- 连续场色标：固定为既有稿件的配色，避免不同环境的可选包改变成图 ----
# 密度、PI 气候态和趋势/差值图分别沿用 viridis、magma、RdBu_r。
SEQ, SEQ_WARM, DIV = "viridis", "magma", "RdBu_r"

# ---- 莫奈风发散色标（仅用于非定量、偏审美强调；玫瑰↔奶白↔青）----
MONET_DIV = LinearSegmentedColormap.from_list(
    "monet_div", ["#7C3E50", "#D98C9C", "#F2ECE3", "#8FBDC4", "#34707F"])

# ---- 小波例外：沿用常规感知顺序色标，不套莫奈 ----
WAVELET_CMAP = "magma"

PROJ = ccrs.PlateCarree()
_GSHHG_LAND = None

# 输出目录（项目根下 figures/），用绝对路径避免 cwd 依赖。
FIGDIR = _LAYOUT_ROOT / 'figures'


def _load_gshhg_land():
    """读取并缓存配置指定的 GSHHG L1 陆地多边形。"""
    global _GSHHG_LAND
    if _GSHHG_LAND is None:
        import geopandas as gpd

        cfg = load_config()
        _GSHHG_LAND = gpd.read_file(cfg["gshhg_path"])
        if _GSHHG_LAND.crs is None:
            _GSHHG_LAND = _GSHHG_LAND.set_crs("EPSG:4326")
        elif _GSHHG_LAND.crs.to_epsg() != 4326:
            _GSHHG_LAND = _GSHHG_LAND.to_crs("EPSG:4326")
    return _GSHHG_LAND





def make_map_ax(ax=None, figsize=(COL1, COL1 * 0.7), *,
                left_labels=True, bottom_labels=True):
    """带 GSHHG 海岸线、经纬网和比例尺的地图 ax（不绘制国境线）。"""
    cfg = load_config()
    r = cfg["regions"]["tc"]   # 地图范围用 TC 统计区
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection=PROJ)
    ax.set_extent([r["lon_min"], r["lon_max"], r["lat_min"], r["lat_max"]], crs=PROJ)
    land = _load_gshhg_land()
    ax.add_geometries(land.geometry, crs=PROJ, facecolor="#EDEAE3",
                      edgecolor="#5A5A5A", linewidth=0.35, zorder=2)
    gl = ax.gridlines(draw_labels=True, lw=0.2, color="#B0B0B0", alpha=0.5)
    gl.top_labels = gl.right_labels = False
    gl.left_labels = left_labels
    gl.bottom_labels = bottom_labels
    return ax


def make_plain_map_ax(ax=None, figsize=(COL1, COL1 * 0.7), *,
                      left_labels=True, bottom_labels=True):
    """Stable Plate Carree map on a standard Matplotlib axis.

    This avoids a Cartopy/Matplotlib multi-GeoAxes colorbar rendering bug while
    retaining the same lon/lat coordinates and extent. Coastlines come from the
    configured GSHHG L1 polygons; administrative borders are intentionally absent.
    """
    cfg = load_config()
    r = cfg["regions"]["tc"]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    xmin, xmax = r["lon_min"], r["lon_max"]
    ymin, ymax = r["lat_min"], r["lat_max"]
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    land = _load_gshhg_land()
    land.cx[xmin:xmax, ymin:ymax].plot(
        ax=ax, facecolor="#EDEAE3", edgecolor="#5A5A5A", linewidth=0.35, zorder=2)
    ax.set_xticks(np.arange(np.ceil(xmin / 20) * 20, xmax + 1, 20))
    ax.set_yticks(np.arange(np.ceil(ymin / 10) * 10, ymax + 1, 10))
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(labelleft=left_labels, labelbottom=bottom_labels)
    ax.grid(color="#B0B0B0", linewidth=0.2, alpha=0.5, zorder=0)
    midlat = 0.5 * (ymin + ymax)
    ax.set_aspect(1.0 / np.cos(np.deg2rad(midlat)), adjustable="box")
    return ax


def panel_letter(ax, s, dx=0.02, dy=0.98):
    """面板字母 a/b/c…（CNS 习惯，左上角加粗）。"""
    ax.text(dx, dy, s, transform=ax.transAxes, fontweight="bold",
            fontsize=13, va="top", ha="left")


def stipple(ax, lon, lat, mask, s=0.5):
    """显著格点打点（趋势 / 合成显著性）。"""
    LON, LAT = np.meshgrid(lon, lat)
    kwargs = {"transform": PROJ} if hasattr(ax, "projection") else {}
    ax.scatter(LON[mask], LAT[mask], s=s, c="#333333", linewidths=0,
               zorder=4, **kwargs)


def add_series(ax, x, y, slope=None, lo=None, hi=None, label=None, color=None):
    """年序列 + Sen 趋势线 + 95% 自助置信带。"""
    color = color or MONET[0]
    ax.plot(x, y, lw=1.4, color=color, label=label)
    if lo is not None and hi is not None:
        ax.fill_between(x, lo, hi, color=color, alpha=0.15, lw=0)
    if slope is not None:
        b = np.mean(y) - slope * np.mean(x)
        xs = np.array([np.min(x), np.max(x)])
        ax.plot(xs, slope * xs + b, "--", lw=1.0, color="#7C3E50")


def save(fig, name):
    """同时输出 600 dpi PNG 与矢量 PDF。"""
    FIGDIR.mkdir(parents=True, exist_ok=True)
    stem = name.rsplit(".", 1)[0]
    # Title-free manuscript variant: captions carry the scientific description,
    # while axes retain labels/units, legends, panel letters and annotations.
    if fig._suptitle is not None:
        fig._suptitle.set_text("")
    for ax in fig.axes:
        for loc in ("left", "center", "right"):
            ax.set_title("", loc=loc)
    # A small, consistent edge prevents clipping while avoiding the large
    # white border introduced by Matplotlib's default tight-layout padding.
    # Keep a restrained outer breathing margin so labels and panel letters do
    # not feel pressed against the page edge in the title-free layout.
    fig.savefig(FIGDIR / f"{stem}.png", dpi=600, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(FIGDIR / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"saved figures/{stem}.png + .pdf")
