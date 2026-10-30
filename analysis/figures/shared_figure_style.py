'Shared agency colours, typography and geographic formatting.'
from __future__ import annotations

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib import font_manager

from matplotlib.ticker import FuncFormatter

import numpy as np

from core.paths import PROJECT,output_path



_names={f.name for f in font_manager.fontManager.ttflist}

FONT=next((f for f in ['Arial','Liberation Sans','DejaVu Sans'] if f in _names),'DejaVu Sans')

plt.rcParams.update({

    'font.family':'sans-serif','font.sans-serif':[FONT,'DejaVu Sans'],

    'mathtext.fontset':'custom','mathtext.rm':FONT,'mathtext.it':FONT+':italic','mathtext.bf':FONT+':bold','mathtext.sf':FONT,

    'font.size':9,'axes.titlesize':9,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,

    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,

    'xtick.major.size':3,'ytick.major.size':3,'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':300})



# Distinct agency colours, fixed across all figures.

COLORS={'PRIMARY':'#222222','USA':'#1F5AA6','TOKYO':'#2E9E44','CMA':'#D1495B','HKO':'#E08E0B'}

LABELS={'PRIMARY':'Primary','USA':'USA','TOKYO':'JMA','CMA':'CMA','HKO':'HKO'}

PERIODS={'early':('#3B6EA8','1966–1995'),'late':('#C8553D','1996–2025')}

LAND='#E9E7E2';COAST='#606060'



def panel(ax,letter,title=''):

    """Bold panel letter followed by a short plain title, left-aligned above the axes."""

    ax.set_title(rf'$\bf{{{letter}}}$   {title}' if title else rf'$\bf{{{letter}}}$',loc='left',pad=4)



_land={}

def land(ax,extent,resolution='i',edge=COAST,lw=.35,face=LAND,zorder=3):

    import geopandas as gpd

    key=resolution

    if key not in _land:

        _land[key]=gpd.read_file(PROJECT/f'data/raw/GSHHG/GSHHS_shp/{resolution}/GSHHS_{resolution}_L1.shp',bbox=(70,-10,190,70))

    _land[key].plot(ax=ax,facecolor=face,edgecolor=edge,linewidth=lw,zorder=zorder)

    ax.set_xlim(extent[0],extent[1]);ax.set_ylim(extent[2],extent[3])



def geo_ticks(ax,xticks,yticks):

    def lon(v,_):

        if np.isclose(abs(v),180):return '180°'

        return f'{abs(v):g}°{"E" if v>0 else "W"}' if not np.isclose(v,0) else '0°'

    def lat(v,_):

        return f'{abs(v):g}°{"N" if v>0 else "S"}' if not np.isclose(v,0) else '0°'

    ax.set_xticks(xticks);ax.set_yticks(yticks)

    ax.xaxis.set_major_formatter(FuncFormatter(lon));ax.yaxis.set_major_formatter(FuncFormatter(lat))

    ax.set_xlabel('');ax.set_ylabel('')

    for s in ['top','right']:ax.spines[s].set_visible(True)



def save(fig,name,number=None,title=None):

    meta={'Creator':'Matplotlib','Author':'Daizhao Huang; Yuanfu Yue'}

    if title:meta['Title']=title

    if number:meta['Subject']=f'Journal of Geophysical Research: Atmospheres Figure {number}'

    fig.savefig(output_path(f'figures/{name}.png'),dpi=300,bbox_inches='tight')

    fig.savefig(output_path(f'figures/{name}.pdf'),bbox_inches='tight',metadata=meta)

    plt.close(fig);print('Figure',name,flush=True)
