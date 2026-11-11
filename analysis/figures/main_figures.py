"""Shared plot data, fonts and geographic axes for the study figures."""
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT


import numpy as np,pandas as pd

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib import font_manager


from matplotlib.text import Text

from matplotlib.ticker import FuncFormatter

from core.paths import PROJECT,output_path

_names={f.name for f in font_manager.fontManager.ttflist}

FONT=next((f for f in ['Arial','Liberation Sans'] if f in _names),'DejaVu Sans')

matplotlib.rcParams.update({'font.family':'sans-serif','font.sans-serif':[FONT,'DejaVu Sans'],'pdf.fonttype':42,'ps.fonttype':42,

    'mathtext.fontset':'custom','mathtext.rm':FONT,'mathtext.it':FONT+':italic','mathtext.bf':FONT+':bold'})

COLORS={'PRIMARY':'#222222','USA':'#1F5AA6','TOKYO':'#2E9E44','CMA':'#D1495B','HKO':'#E08E0B'}

DISPLAY={'PRIMARY':'Primary','USA':'USA','TOKYO':'JMA','CMA':'CMA','HKO':'HKO'}

AG3=['USA','TOKYO','CMA'];AG4=AG3+['HKO']


ZA=_LAYOUT_ROOT / 'data/reference/agency_comparison';ZV=_LAYOUT_ROOT / 'data/reference/environment_validation'

H=lambda n:pd.read_csv(output_path('results/agency_comparison/figure_series/'+n))

def scale_typography(fig,scale):

    """Scale figure lettering with the available Arial-compatible font."""

    fig.canvas.draw()

    for text in fig.findobj(match=Text):

        if text.get_visible() and text.get_fontsize()>0:

            text.set_fontfamily(FONT);text.set_fontsize(text.get_fontsize()*scale)

    fig.canvas.draw()

def add_coast(ax,extent):

    import geopandas as gpd

    land=gpd.read_file(PROJECT/'data/raw/GSHHG/GSHHS_shp/i/GSHHS_i_L1.shp',bbox=extent)

    land.plot(ax=ax,facecolor='#E6E4DF',edgecolor='#505050',linewidth=0.35,zorder=3)

    ax.set_xlim(extent[0],extent[2]);ax.set_ylim(extent[1],extent[3])

    # Map scale true at the central latitude of the panel (same convention as Figure 1).

    ax.set_aspect(1/np.cos(np.deg2rad((extent[1]+extent[3])/2)),adjustable='box')

def format_geo_axes(ax):

    lon=lambda v,_:'180°' if np.isclose(abs(v),180) else ('0°' if np.isclose(v,0) else f"{abs(v):g}°{'E' if v>0 else 'W'}")

    lat=lambda v,_:'0°' if np.isclose(v,0) else f"{abs(v):g}°{'N' if v>0 else 'S'}"

    ax.xaxis.set_major_formatter(FuncFormatter(lon));ax.yaxis.set_major_formatter(FuncFormatter(lat));ax.set_xlabel('');ax.set_ylabel('')

def open_panel_axis(ax,*,right=False):

    ax.spines['top'].set_visible(False)

    if right:

        ax.spines['left'].set_visible(False);ax.yaxis.set_label_position('right');ax.yaxis.tick_right()

        ax.tick_params(axis='y',left=False,labelleft=False,right=True,labelright=True)

    else:ax.spines['right'].set_visible(False)
