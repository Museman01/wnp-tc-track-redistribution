'Si four agency figures.'
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np,pandas as pd

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

from matplotlib.patches import Patch

from core.paths import PROJECT,output_path

from analysis.figures.main_figures import COLORS, DISPLAY, AG4, scale_typography, ZA



R=lambda n:pd.read_csv(output_path('results/agency_comparison/'+n))

OUTPUT=None
NUM={'S01':1,'S02':2,'S03':3,'S06':6,'S07':7,'S09':9,'S11':11,'S12':12,'S13':13}



def letters(axes,x=0.015,y=0.985,size=12):

    for label,ax in zip('abcdefghijklmnopqrstuvwxyz',np.ravel(axes)):

        ax.text(x,y,label,transform=ax.transAxes,va='top',ha='left',fontsize=size,fontweight='bold',zorder=10)



def open_axis(ax,right=False):

    ax.spines['top'].set_visible(False)

    if right:

        ax.spines['left'].set_visible(False);ax.yaxis.set_label_position('right');ax.yaxis.tick_right()

        ax.tick_params(axis='y',left=False,labelleft=False,right=True,labelright=True)

    else:ax.spines['right'].set_visible(False)



def save(fig,stem,key,title,scale):

    if OUTPUT is None: raise ValueError('Provide an output directory')
    if any((OUTPUT / (stem+ext)).exists() for ext in ('.png','.pdf')): raise FileExistsError('Refusing to overwrite an existing figure')
    scale_typography(fig,scale)

    meta={'Title':title,'Subject':f'Journal of Geophysical Research: Atmospheres Supporting Information Figure S{NUM[key]}','Creator':'Matplotlib'}

    fig.savefig(OUTPUT / f'{stem}.png',dpi=300,bbox_inches='tight')

    fig.savefig(OUTPUT / f'{stem}.pdf',bbox_inches='tight',metadata=meta);plt.close(fig);print('Figure',stem,flush=True)



def agency_legend(ax,**kw):

    ax.legend(handles=[Line2D([0],[0],color=COLORS[a],lw=1.5,label=DISPLAY[a]) for a in AG4],frameon=False,**kw)



# ------------------------------------------------------------------ S1

def S01():

    c=R('coverage/coverage_four_agency.csv')

    fig,axes=plt.subplots(3,1,figsize=(8.7,8.0),sharex=True,constrained_layout=True,gridspec_kw={'height_ratios':[1,1,1]})

    for a in AG4:

        g=c[c.agency.eq(a)]

        axes[0].plot(g.season,g.n_storms,color=COLORS[a],lw=1.2)

        axes[1].plot(g.season,g.n_positions,color=COLORS[a],lw=1.2)

        axes[2].plot(g.season,g.wind_available_pct,color=COLORS[a],lw=1.2)

    axes[0].set_ylabel('Storms meeting native\nTS criterion');axes[1].set_ylabel('Native six-hourly\npositions')

    axes[2].set(ylabel='Wind-field availability (%)',xlabel='Year');axes[2].set_ylim(-5,105)

    agency_legend(axes[0],ncol=4,loc='upper right')

    for ax in axes:open_axis(ax)

    letters(axes)

    save(fig,'FigS01_agency_coverage_1945_2025','S01','Native-record coverage',1.10)



# ------------------------------------------------------------------ S2

def S02():

    p=R('coverage/projection_1945_2025_four_agency.csv')

    fig,ax=plt.subplots(figsize=(8.2,6.8),constrained_layout=True)

    for a in AG4:

        g=p[p.agency.eq(a)];ax.plot(g.year,g.projection,color=COLORS[a],lw=1.2)

    ax.axhline(0,color='0.45',lw=0.8);ax.axvline(1966,color='0.45',ls=':',lw=1.0);ax.axvline(1996,color='0.45',ls='--',lw=1.0)

    ax.set(xlabel='Year',ylabel='Projection onto the 1966–2025 contrast')

    agency_legend(ax,loc='upper left');open_axis(ax)

    save(fig,'FigS02_early_year_projection','S02','Early-year projection',1.10)



# ------------------------------------------------------------------ S3

def S03():

    s=R('coverage/start_year_four_agency.csv')

    fig,axes=plt.subplots(1,2,figsize=(11,5.4),constrained_layout=True)

    starts=[1945,1966,1982];w=0.2

    for i,a in enumerate(AG4):

        g=s[s.agency.eq(a)].set_index('start')

        for ax,col in [(axes[0],'correlation_with_1966'),(axes[1],'total_variation')]:

            x=[k+(i-1.5)*w for k,st in enumerate(starts) if st in g.index];v=[g.loc[st,col] for st in starts if st in g.index]

            ax.bar(x,v,width=w,color=COLORS[a])

    axes[0].set(xticks=range(3),xticklabels=starts,xlabel='Start year',ylabel='Spatial correlation with 1966 definition')

    axes[1].set(xticks=range(3),xticklabels=starts,xlabel='Start year',ylabel='Total variation distance')

    open_axis(axes[0]);open_axis(axes[1],right=True)

    fig.legend(handles=[Patch(color=COLORS[a],label=DISPLAY[a]) for a in AG4],frameon=False,ncol=4,loc='upper center',bbox_to_anchor=(0.5,1.07))

    letters(axes)

    save(fig,'FigS03_start_year_sensitivity','S03','Start-year sensitivity',1.20)



# ------------------------------------------------------------------ S6

def S06():

    z=pd.read_csv(ZA/'04_track_density_sensitivity/path_definition_sensitivity.csv').query('grid_deg == 2.5 and block_years == 3')

    h=R('path_definition/path_definition_hko.csv').query('grid_deg == 2.5 and block_years == 3')

    path=pd.concat([z[z.agency.isin(['USA','TOKYO','CMA'])],h])

    defs=['track_point','storm_equal','binary_occupancy','line_length','pre_landfall_track_point','full_life_track_point']

    labels=['Points','Storm equal','Binary','Line length','Pre-landfall','Full life']

    fig,axes=plt.subplots(1,2,figsize=(11,4.3),constrained_layout=True);w=0.2;base=np.arange(len(defs))

    for i,a in enumerate(AG4):

        part=path[path.agency.eq(a)].set_index('definition').loc[defs]

        axes[0].bar(base+(i-1.5)*w,part.total_variation,width=w,color=COLORS[a],label=DISPLAY[a])

        axes[1].bar(base+(i-1.5)*w,part.block_permutation_p,width=w,color=COLORS[a])

    axes[0].set(xticks=base,xticklabels=labels,ylabel='Total variation distance')

    axes[1].set(xticks=base,xticklabels=labels,ylabel='Block-permutation p');axes[1].axhline(.05,color='0.45',ls=':',lw=1);axes[1].set_ylim(0,.055)

    for ax in axes:ax.tick_params(axis='x',rotation=25)

    open_axis(axes[0]);open_axis(axes[1],right=True)

    fig.legend(handles=[Patch(color=COLORS[a],label=DISPLAY[a]) for a in AG4],frameon=False,ncol=4,loc='upper center',bbox_to_anchor=(0.5,1.07))

    letters(axes)

    save(fig,'FigS06_path_definition_sensitivity','S06','Path-definition sensitivity',1.22)



# ------------------------------------------------------------------ S7

def S07():

    c=R('climate/climate_mode_four_agency.csv').set_index('agency').loc[AG4]

    fig,axes=plt.subplots(1,2,figsize=(8.5,3.8),constrained_layout=True);x=np.arange(4)

    axes[0].bar(x,c.partial_r2,color=[COLORS[a] for a in AG4]);axes[1].bar(x,c.block_permutation_p,color=[COLORS[a] for a in AG4])

    axes[0].set(xticks=x,xticklabels=[DISPLAY[a] for a in AG4],ylabel='Partial R²')

    axes[1].set(xticks=x,xticklabels=[DISPLAY[a] for a in AG4],ylabel='Block-permutation p');axes[1].axhline(.05,color='0.45',ls=':',lw=1)
    axes[1].set_yscale('log');axes[1].set_ylim(1e-4,0.1)   # p values of 0.0005-0.0033 stay visible against the 0.05 line
    axes[0].set_ylim(top=axes[0].get_ylim()[1]*1.12)   # headroom for the panel letter

    open_axis(axes[0]);open_axis(axes[1],right=True);letters(axes)

    save(fig,'FigS07_climate_mode_adjustment','S07','Climate-mode adjustment',1.20)



# ------------------------------------------------------------------ S9

def S09():

    an=R('speed/speed_decomposition_annual_four_agency.csv');tr=R('speed/speed_decomposition_trends_four_agency.csv')

    fig,axes=plt.subplots(2,2,figsize=(11,7.6),constrained_layout=True)

    tot=an[an.speed_type.eq('total')]

    for a in AG4:

        g=tot[tot.agency.eq(a)];axes[0,0].plot(g.season,g.basin_speed,color=COLORS[a],lw=1.2,label=DISPLAY[a])

    axes[0,0].set(ylabel='Basin-mean speed (km h$^{-1}$)');agency_legend(axes[0,0],ncol=4,loc='lower center',bbox_to_anchor=(0.5,1.0))

    terms=[('track_shift','Track-density term','#4C78A8'),('regional_speed','Regional-speed term','#E45756'),('covariance','Covariance term','#72B7B2')]

    j=tot[tot.agency.eq('TOKYO')]

    for v,l,c in terms:axes[0,1].plot(j.season,j[v],color=c,lw=1.2,label=l)

    axes[0,1].axhline(0,color='0.45',lw=0.8);axes[0,1].set(ylabel='Speed anomaly (km h$^{-1}$)');axes[0,1].legend(frameon=False,ncol=3,loc='lower center',bbox_to_anchor=(0.5,1.0),handlelength=1.6,columnspacing=1.0)   # also the key for panel c

    comp=tr[tr.speed_type.eq('total')&tr.dataset.isin([f'{a}_feng' for a in AG4])]

    w=0.24;x=np.arange(4)

    for k,(v,l,c) in enumerate(terms):

        part=comp[comp.variable.eq(v)].set_index('dataset').loc[[f'{a}_feng' for a in AG4]]

        vals=part.period_change_percent_of_mean_basin_speed.to_numpy();axes[1,0].bar(x+(k-1)*w,vals,width=w,color=c,label=l)

        for i,(val,q) in enumerate(zip(vals,part.mk_p_fdr_bh)):

            if q<0.05:axes[1,0].text(x[i]+(k-1)*w,val+0.15,'*',ha='center',va='bottom',fontsize=12)

    axes[1,0].axhline(0,color='0.4',lw=0.8);axes[1,0].set(xticks=x,xticklabels=[DISPLAY[a] for a in AG4],ylabel='Change (% of mean basin speed)')

    axes[1,0].set_ylim(top=9.0)

    pr=an[an.agency.eq('USA_project')]

    ex=pd.read_csv(PROJECT/'data/processed/annual_analysis_table.csv');ex=ex[ex.scope.eq('annual')&ex.season.between(1966,2025)]

    axes[1,1].plot(pr.season,pr.basin_speed,color='#1F5AA6',lw=1.2,label='Point-weighted decomposition')

    axes[1,1].plot(ex.season,ex.speed,color='#222222',lw=1.2,label='Storm-equal metric')

    axes[1,1].set(xlabel='Year',ylabel='Basin-mean speed (km h$^{-1}$)');axes[1,1].set_ylim(top=axes[1,1].get_ylim()[1]+0.30*np.ptp(axes[1,1].get_ylim()));axes[1,1].legend(frameon=False,ncol=1,loc='upper right')

    for ax in axes[:,0]:open_axis(ax)

    for ax in axes[:,1]:open_axis(ax,right=True)

    letters(axes)

    save(fig,'FigS09_translation_speed_decomposition','S09','Translation-speed decomposition',1.15)



# ------------------------------------------------------------------ S11

def S11():

    ann=R('landfall/landfall_annual_four_agency.csv');s=R('landfall/landfall_summary_four_agency.csv')

    s=s[s.metric.eq('mean_lat')&s.start.eq(1966)]

    fig,axes=plt.subplots(2,2,figsize=(11.4,8.0),constrained_layout=True)

    for ax,d in [(axes[0,0],'first_landfall'),(axes[0,1],'all_events')]:

        for a in AG4:

            g=ann[ann.agency.eq(a)&ann.definition.eq(d)];ax.plot(g.year,g.mean_lat,color=COLORS[a],lw=1.1)

        ax.set_ylabel('Annual mean landfall latitude (°N)')

    agency_legend(axes[0,0],ncol=4,loc='lower center')

    e=s[s.end.eq(2025)];rows=[(a,d) for d in ['first_landfall','all_events'] for a in AG4]

    for i,(a,d) in enumerate(rows):

        r=e[e.agency.eq(a)&e.definition.eq(d)].iloc[0]

        axes[1,0].errorbar(i,r.period_difference,yerr=[[r.period_difference-r.period_ci_low],[r.period_ci_high-r.period_difference]],fmt='o',color='0.2',ms=6,capsize=4)

    axes[1,0].axhline(0,color='0.5',lw=0.8)

    axes[1,0].set(xticks=range(8),xticklabels=[f"{DISPLAY[a]}\n{'First' if d=='first_landfall' else 'All'}" for a,d in rows],ylabel='Later-minus-earlier latitude (°)')

    f=s[s.definition.eq('first_landfall')]

    for a in AG4:

        g=f[f.agency.eq(a)].sort_values('end');axes[1,1].plot(g.end,g.period_difference,'-o',color=COLORS[a],mfc=COLORS[a],ms=6,label=DISPLAY[a])

    axes[1,1].set(xticks=[2024,2025],xlim=(2023.8,2025.2),xlabel='Analysis end year',ylabel='First-landfall latitude change (°)')

    axes[1,1].legend(frameon=False,ncol=4,loc='upper right')

    for ax in axes[:,0]:open_axis(ax)

    for ax in axes[:,1]:open_axis(ax,right=True)

    for ax in axes[0]:ax.set_ylim(bottom=ax.get_ylim()[0]-1.0)

    letters(axes)

    save(fig,'FigS11_landfall_latitude_diagnostic','S11','Landfall-latitude diagnostics',1.15)



# ------------------------------------------------------------------ S12

def S12():

    c=R('sensitivity/cutpoint_four_agency.csv');c=c[c.scheme.eq('full_record')]

    mets=[('path_total_variation','Total variation distance'),('first_landfall_latitude_difference','Latitude change (°)'),('north_share_difference','Share change (percentage points)'),('lmi_latitude_difference','Latitude change (°)')]

    fig,axes=plt.subplots(4,2,figsize=(9.2,13.2),sharex=True,constrained_layout=True)

    for r,(m,lab) in enumerate(mets):

        for a in AG4:

            g=c[c.agency.eq(a)&c.metric.eq(m)].sort_values('cutpoint')

            axes[r,0].plot(g.cutpoint,g.effect*(100 if m=='north_share_difference' else 1),'-o',color=COLORS[a],ms=3.5,lw=1.2,label=DISPLAY[a])

            axes[r,1].plot(g.cutpoint,g.p_value,'-o',color=COLORS[a],ms=3.5,lw=1.2)

        axes[r,0].set_ylabel(lab);axes[r,1].set_ylabel('Block-permutation p');axes[r,1].axhline(0.05,color='0.45',ls=':',lw=1)

        for ax in axes[r]:ax.axvline(1996,color='0.4',ls='--',lw=0.9)

        for ax in axes[r]:lo_,hi_=ax.get_ylim();ax.set_ylim(lo_,hi_+0.12*(hi_-lo_))   # headroom for the panel letters
        open_axis(axes[r,0]);open_axis(axes[r,1],right=True)

    axes[0,0].legend(frameon=False,loc='upper left',bbox_to_anchor=(0.08,1.0))


    for ax in axes[-1]:ax.set_xlabel('First year of late period')

    letters(np.r_[axes[:,0],axes[:,1]])

    save(fig,'FigS12_cutpoint_sensitivity','S12','Moving-cutpoint sensitivity',1.40)



# ------------------------------------------------------------------ S13

def S13():

    t=R('sensitivity/coast_grouping_four_agency.csv').query("taiwan_rule == 'taiwan_north' and denominator == 'named_only'")

    fig,axes=plt.subplots(3,1,figsize=(8.6,9.8),sharex=True,constrained_layout=True)

    for ax,rule in zip(axes,['all_events','first_any','first_named']):

        for a in AG4:

            g=t[t.event_rule.eq(rule)&t.agency.eq(a)].sort_values('threshold_km');sig=g.q_bh<0.05

            ax.plot(g.threshold_km,g.change_percentage_points,color=COLORS[a],lw=1.5)

            ax.scatter(g.threshold_km[sig],g.change_percentage_points[sig],s=42,facecolors=COLORS[a],edgecolors=COLORS[a],zorder=3)

            ax.scatter(g.threshold_km[~sig],g.change_percentage_points[~sig],s=42,facecolors='white',edgecolors=COLORS[a],zorder=3)

        ax.set_ylabel('North-share change\n(percentage points)');open_axis(ax)

        lo,hi=ax.get_ylim();ax.set_ylim(lo-0.4,hi+0.4)

    axes[-1].set(xticks=[25,50,75,100],xlabel='Coastline-segment threshold (km)')

    agency_legend(axes[0],ncol=4,loc='upper right');letters(axes)

    save(fig,'FigS13_exclusive_coast_threshold_sensitivity','S13','Coast-threshold sensitivity',1.12)



if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('figures',nargs='*',choices=list(NUM))
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    if not args.execute: parser.error('add --execute to create figures')
    OUTPUT=(PROJECT / args.output_dir).resolve()
    if not OUTPUT.is_relative_to(PROJECT.resolve()): parser.error('Output directory must remain within the project')
    if OUTPUT.exists(): parser.error('Provide a new output directory')
    OUTPUT.mkdir(parents=True)
    for key in (args.figures or list(NUM)): globals()[key]()
