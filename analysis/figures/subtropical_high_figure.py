"""Redraw Figure S9 at 70% height from the existing annual series and trend table."""
from __future__ import annotations
from analysis.figures.figure_raster import render_png, apply_tnr
from analysis.figures.si_page_type import fit_si_type
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core.paths import PROJECT
import analysis.figures.main_figures as mf


def render(stem):
    stem=(PROJECT/Path(stem)).resolve()
    if not stem.is_relative_to(PROJECT.resolve()):raise ValueError('Output must remain within the project')
    if any(Path(str(stem)+s).exists() for s in ('.png','.pdf')):raise FileExistsError(stem)
    annual=pd.read_csv(PROJECT/'data/processed/subtropical_high_indices.csv').query('1966 <= season <= 2025')
    trends=pd.read_csv(PROJECT/'data/reference/environment_validation/circulation_metric_validation/results/wnpsh_metric_trends.csv')
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':[mf.FONT,'DejaVu Sans'],'pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(694.7995/72,349.9148/72),sharex=True)
    fig.subplots_adjust(left=.1077,right=.9193,bottom=.038,top=.845,wspace=.08,hspace=.34)
    cols=[('wpsh_area','WNPSH area\n($10^6$ km²)',np.arange(0,18,5)),
          ('wpsh_intensity','WNPSH\nintensity',np.arange(0,41,10)),
          ('ridge_line','Ridge latitude (°N)',np.arange(21,28,2)),
          ('west_ridge_point','Western ridge\npoint (°E)',np.arange(110,161,20))]
    years=np.arange(1970,2030,10)
    for i,(ax,(column,label,ticks)) in enumerate(zip(axes.flat,cols)):
        old=annual[annual.season<1982]; modern=annual[annual.season>=1982]
        ax.plot(old.season,old[column],color='#adadad',lw=1.15,zorder=2)
        ax.plot(modern.season,modern[column],color='#4c72b0',lw=1.45,zorder=3)
        ax.axvline(1982,color='#777777',linestyle=':',lw=.9,zorder=1)
        valid=modern[['season',column]].dropna()
        rows=trends[(trends.variable==column)&(trends.start==1982)&(trends.end==2025)]
        if len(rows)!=1 or len(valid)!=int(rows.iloc[0]['n']):raise ValueError('Trend input mismatch: '+column)
        slope=float(rows.iloc[0].sen_slope_per_year)
        x=np.array([valid.season.min(),valid.season.max()],float)
        y=valid[column].mean()+slope*(x-valid.season.mean())
        ax.plot(x,y,color='#7c3e50',ls='--',lw=1.15,zorder=4)
        ax.set_xlim(1960,2025);ax.set_xticks(years);ax.set_yticks(ticks)
        ax.set_ylabel(label,fontsize=12.65,labelpad=4)
        ax.tick_params(axis='both',labelsize=12.65,width=.8,length=3.5,pad=2)
        if i % 2:
            low, high = ax.get_ylim()
            ax.set_ylim(low, high + 0.18 * (high - low))
        ax.text(.985 if i % 2 else .015, .985, 'abcd'[i], transform=ax.transAxes,
                va='top', ha='right' if i % 2 else 'left', fontsize=12.65,
                fontweight='bold', zorder=8)
        ax.spines['top'].set_visible(i<2)
        ax.spines['bottom'].set_visible(i>=2)
        if i%2:
            ax.spines['left'].set_visible(False);ax.spines['right'].set_visible(True)
            ax.yaxis.tick_right();ax.yaxis.set_label_position('right')
        else:
            ax.spines['left'].set_visible(True);ax.spines['right'].set_visible(False)
        if i<2:
            ax.xaxis.set_label_position('top');ax.set_xlabel('Year',fontsize=12.65,labelpad=8)
            ax.tick_params(axis='x',top=True,labeltop=True,bottom=False,labelbottom=False)
        else:
            ax.tick_params(axis='x',top=False,labeltop=False,bottom=True,labelbottom=False)
    stem.parent.mkdir(parents=True,exist_ok=True)
    fit_si_type(fig)
    from matplotlib.transforms import blended_transform_factory
    fig.canvas.draw()
    for column in range(2):
        upper, lower = axes[0, column], axes[1, column]
        middle = (upper.get_position().y0 + lower.get_position().y1) / 2
        upper.text(1982, middle, '1982', transform=blended_transform_factory(upper.transData, fig.transFigure), ha='center', va='center', fontsize=16, clip_on=False)
    apply_tnr(fig)
    fig.savefig(str(stem)+'.pdf',bbox_inches='tight')
    render_png(str(stem)+'.pdf')
    plt.close(fig)
    print('Rendered',stem)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args();render(args.output)
