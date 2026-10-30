'Main figures for tropical cyclone track and landfall analyses.'
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import sys,shutil

import numpy as np,pandas as pd

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib import font_manager

from matplotlib.lines import Line2D

from matplotlib.text import Text

from matplotlib.ticker import FuncFormatter

from matplotlib.patches import FancyBboxPatch

from matplotlib.colors import TwoSlopeNorm

from core.paths import PROJECT,output_path



_names={f.name for f in font_manager.fontManager.ttflist}

FONT=next((f for f in ['Arial','Liberation Sans'] if f in _names),'DejaVu Sans')

matplotlib.rcParams.update({'font.family':'sans-serif','font.sans-serif':[FONT,'DejaVu Sans'],'pdf.fonttype':42,'ps.fonttype':42,

    'mathtext.fontset':'custom','mathtext.rm':FONT,'mathtext.it':FONT+':italic','mathtext.bf':FONT+':bold'})





COLORS={'PRIMARY':'#222222','USA':'#1F5AA6','TOKYO':'#2E9E44','CMA':'#D1495B','HKO':'#E08E0B'}

DISPLAY={'PRIMARY':'Primary','USA':'USA','TOKYO':'JMA','CMA':'CMA','HKO':'HKO'}

AG3=['USA','TOKYO','CMA'];AG4=AG3+['HKO']

Z=_LAYOUT_ROOT / 'data/reference'

ZA=_LAYOUT_ROOT / 'data/reference/agency_comparison';ZV=_LAYOUT_ROOT / 'data/reference/environment_validation'

H=lambda n:pd.read_csv(output_path('results/agency_comparison/figure_series/'+n))



def scale_typography(fig,scale):

    """As figure_typography.scale_figure_typography, with the available Arial-compatible font."""

    fig.canvas.draw()

    for text in fig.findobj(match=Text):

        if text.get_visible() and text.get_fontsize()>0:

            text.set_fontfamily(FONT);text.set_fontsize(text.get_fontsize()*scale)

    fig.canvas.draw()



def save(fig,stem,*,scale,dpi,number,title):

    scale_typography(fig,scale)

    meta={'Title':title,'Subject':f'Journal of Geophysical Research: Atmospheres Figure {number}','Creator':'Matplotlib'}

    fig.savefig(output_path(f'figures/{stem}.png'),dpi=dpi,bbox_inches='tight')

    fig.savefig(output_path(f'figures/{stem}.pdf'),bbox_inches='tight',metadata=meta);plt.close(fig);print('Figure',stem,flush=True)



def panel_labels(axes,size=11):

    for label,ax in zip('abcdefghijklmnopqrstuvwxyz',np.asarray(axes).flat):

        ax.text(0.01,0.98,label,transform=ax.transAxes,ha='left',va='top',fontsize=size,fontweight='bold',zorder=10)



def open_panel_axes(axes,top_row_on_top=True):

    """Panel layout: right column carries its y axis on the right; top-row x axes at the top."""

    for r,row in enumerate(np.atleast_2d(axes)):

        for column,ax in enumerate(row):

            if r==0 and top_row_on_top:

                ax.xaxis.tick_top();ax.xaxis.set_label_position('top');ax.spines['bottom'].set_visible(False)

            else:ax.spines['top'].set_visible(False)

            if column==0:ax.spines['right'].set_visible(False)

            else:

                ax.spines['left'].set_visible(False);ax.spines['right'].set_visible(True)

                ax.yaxis.set_label_position('right');ax.yaxis.tick_right()

                ax.tick_params(axis='y',left=False,labelleft=False,right=True,labelright=True)



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



def panel_letters39(axes):

    for label,ax in zip('abcdefghijklmnopqrstuvwxyz',np.ravel(axes)):

        ax.text(0.015,0.985,label,transform=ax.transAxes,va='top',ha='left',fontsize=12,fontweight='bold',zorder=10)



def open_panel_axis(ax,*,right=False):

    ax.spines['top'].set_visible(False)

    if right:

        ax.spines['left'].set_visible(False);ax.yaxis.set_label_position('right');ax.yaxis.tick_right()

        ax.tick_params(axis='y',left=False,labelleft=False,right=True,labelright=True)

    else:ax.spines['right'].set_visible(False)



# ------------------------------------------------------------------ Figure 2

def figure2():

    annual=pd.read_csv(ZA/'03_common_storms/core_crossagency_annual.csv')

    red=pd.read_csv(ZA/'03_common_storms/core_crossagency_recheck_redistribution.csv')

    land=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/landfall_latitude_summary.csv')

    land=land.loc[land.start.eq(1966)&land.end.eq(2025)&land.metric.eq('mean_lat')&land.definition.eq('first_landfall')]

    hf=H('frequency.csv');hl=H('lmi_common.csv');tv=pd.read_csv(_LAYOUT_ROOT / 'results/tracks/hko_track_field_tests.csv').query("sample=='native_all' and agency=='HKO'").iloc[0]

    hland=H('landfall.csv').query("agency=='HKO' and definition=='first_landfall'").iloc[0]

    fig,axes=plt.subplots(2,2,figsize=(11.4,7.6),constrained_layout=True)

    for agency in AG3:

        part=annual.loc[annual.agency.eq(agency)]

        axes[0,0].plot(part.season,part.n_tc,color=COLORS[agency],lw=1.3,label=DISPLAY[agency])

    axes[0,0].plot(hf.season,hf.n_tc,color=COLORS['HKO'],lw=1.3,label='HKO')

    # Panel b: storms common to all four records (1,465 storms).

    l4=pd.read_csv(output_path('results/agency_comparison/lmi_common_four_agency_annual.csv'))

    for agency in AG4:

        part=l4.loc[l4.agency.eq(agency)];axes[0,1].plot(part.season,part.mean_lmi_lat_common,color=COLORS[agency],lw=1.3,label=DISPLAY[agency])

    axes[0,0].set(xlabel='Year',ylabel='TC count');axes[0,1].set(xlabel='Year',ylabel='Mean LMI\nlatitude (°N)')

    axes[0,0].legend(frameon=False,ncol=4)

    x=np.arange(4);r=red.set_index('agency').loc[AG3]

    tvs=list(r.total_variation)+[tv.total_variation];ps=list(r.block_permutation_p)+[tv.block_permutation_p]

    bars=axes[1,0].bar(x,tvs,color=[COLORS[a] for a in AG4])

    axes[1,0].set(xticks=x,xticklabels=[DISPLAY[a] for a in AG4],ylabel='Total-variation\ndistance')

    for bar,p in zip(bars,ps):axes[1,0].text(bar.get_x()+bar.get_width()/2,bar.get_height()+.003,f'p={p:.4f}',ha='center',fontsize=8)

    l3=land.set_index('agency').loc[AG3]

    eff=np.r_[l3.period_difference.values,hland.period_difference];lo=np.r_[l3.period_ci_low.values,hland.period_ci_low]

    hi=np.r_[l3.period_ci_high.values,hland.period_ci_high];p=np.r_[l3.period_block_p.values,hland.period_block_p]

    axes[1,1].errorbar(x,eff,yerr=[eff-lo,hi-eff],fmt='none',ecolor='0.25',capsize=4)

    axes[1,1].scatter(x,eff,s=55,c=[COLORS[a] for a in AG4],zorder=3);axes[1,1].axhline(0,color='0.45',lw=.8)

    axes[1,1].set(xticks=x,xticklabels=[DISPLAY[a] for a in AG4],ylabel='First-landfall latitude\nchange (°)')

    axes[1,1].set_xlim(-0.3,3.3);axes[1,1].set_ylim(-0.1,float(hi.max()+0.35))

    for i in range(4):axes[1,1].text(i,hi[i]+.08,f'p={p[i]:.4f}',ha='center',fontsize=7.5)

    open_panel_axes(axes);panel_labels(axes)

    save(fig,'Fig02_cross_agency_core',scale=1.32,dpi=300,number=2,title='Cross-agency core results')



# ------------------------------------------------------------------ Figure 3

def figure3():

    f=np.load(_LAYOUT_ROOT / 'results/genesis/PRIMARY_fields.npz')

    annual=pd.read_csv(ZV/'supplemental_validation/results/wnp_tc_redistribution_index_annual.csv')

    robust=pd.read_csv(ZV/'supplemental_validation/results/wnp_tc_robustness_matrix.csv')

    lon_edges=f['lon_edges'];lat_edges=f['lat_edges'];change=f['total_field'].reshape(len(lat_edges)-1,len(lon_edges)-1)*100

    vmax=np.nanmax(np.abs(change))

    fig=plt.figure(figsize=(12.2,7.3),constrained_layout=True)

    gs=fig.add_gridspec(2,3,height_ratios=[0.86,1.0],width_ratios=[0.055,1.0,1.0])

    cax=fig.add_subplot(gs[0,0]);axa=fig.add_subplot(gs[0,1]);axb=fig.add_subplot(gs[0,2]);axc=fig.add_subplot(gs[1,:])

    m=axa.pcolormesh(lon_edges,lat_edges,change,cmap='RdBu_r',vmin=-vmax,vmax=vmax,shading='flat')

    add_coast(axa,(100,0,180,40));format_geo_axes(axa)

    cb=fig.colorbar(m,cax=cax,orientation='vertical');cb.set_label('1996–2025 minus 1966–1995\nshare (percentage\npoints per cell)',labelpad=9)

    cax.yaxis.set_ticks_position('left');cax.yaxis.set_label_position('left')

    for agency in ['PRIMARY']+AG3:

        d=annual[(annual.agency==agency)&(annual.weighting=='track_point')].sort_values('year')

        axb.plot(d.year,d.index_oos,lw=1.5 if agency=='PRIMARY' else 1.1,color=COLORS[agency],label=DISPLAY[agency],zorder=4 if agency=='PRIMARY' else 3)

    h=H('index.csv');axb.plot(h.year,h.index_oos,lw=1.1,color=COLORS['HKO'],label='HKO')

    axb.axhline(0,color='0.45',lw=0.7);axb.axvline(1995.5,color='0.55',ls='--',lw=0.8)

    axb.set(xlabel='Year',ylabel='Leave-one-year-out\nredistribution index');axb.legend(frameon=False,ncol=2,loc='upper left',bbox_to_anchor=(0.04,1.0))

    axb.set_box_aspect(40/80/np.cos(np.deg2rad(20)));axa.set_xticks([100,120,140,160,180]);axa.set_yticks([0,10,20,30,40])

    selected=[]

    for label,query in [

        ('1° grid',"analysis == 'grid_primary' and catalog == 'PRIMARY' and weighting == 'track_point' and grid_deg == 1.0"),

        ('2.5° grid\n(main)',"analysis == 'grid_primary' and catalog == 'PRIMARY' and weighting == 'track_point' and grid_deg == 2.5"),

        ('5° grid',"analysis == 'grid_primary' and catalog == 'PRIMARY' and weighting == 'track_point' and grid_deg == 5.0"),

        ('Later period\nends 2024',"analysis == 'end_2024_drop_1995' and catalog == 'PRIMARY' and weighting == 'track_point'"),

        ('Earlier period\nstarts 1967',"analysis == 'start_1967_drop_1996' and catalog == 'PRIMARY' and weighting == 'track_point'"),

        ('2020–2025\nexcluded',"analysis == 'exclude_2020_2025' and catalog == 'PRIMARY' and weighting == 'track_point'"),

        ('Typhoons only,\ntrack-point\nweighting',"analysis == 'typhoon_threshold' and catalog == 'PRIMARY_TY' and weighting == 'track_point'"),

        ('Typhoons only,\nequal storm\nweighting',"analysis == 'typhoon_threshold' and catalog == 'PRIMARY_TY' and weighting == 'storm_equal'")]:

        selected.append((label,robust.query(query).iloc[0].block_permutation_p))

    labels,vals=zip(*selected)

    x=np.arange(len(labels));bars=axc.bar(x,vals,color=['#8DA0CB' if p>=0.05 else '#4C72B0' for p in vals])

    axc.axhline(0.05,color='#B44A4A',ls='--',lw=1.0)

    axc.set_xticks(x,labels);axc.set_ylabel('Block-permutation p')

    from matplotlib.patches import Patch

    axc.legend(handles=[Patch(color='#4C72B0',label='p < 0.05'),Patch(color='#8DA0CB',label='p ≥ 0.05'),Line2D([0],[0],color='#B44A4A',ls='--',lw=1.0,label='p = 0.05')],frameon=False,loc='upper left',bbox_to_anchor=(0.03,1.0))

    for bar,p in zip(bars,vals):axc.text(bar.get_x()+bar.get_width()/2,p+0.006,f'{p:.3f}',ha='center',fontsize=8)

    open_panel_axis(axb,right=True);open_panel_axis(axc)

    # Colour bar spans the map height (layout frozen after the constrained pass).

    panel_letters39([axa,axb,axc]);scale_typography(fig,1.20);fig.canvas.draw();fig.set_layout_engine('none');pa=axa.get_position();pc=cax.get_position();cax.set_position([pc.x0,pa.y0,pc.width,pa.height])

    save(fig,'Fig03_track_redistribution_index_robustness',scale=1.0,dpi=320,number=3,title='Track redistribution, annual index and robustness')



# ------------------------------------------------------------------ Figure 4

def figure4():

    f=np.load(_LAYOUT_ROOT / 'results/genesis/PRIMARY_fields.npz')

    cats=['PRIMARY','USA','TOKYO','CMA','HKO']

    s=pd.concat([pd.read_csv(output_path(f'results/genesis/{c}_summary.csv')) for c in cats]).set_index('catalog').loc[cats]

    lon_edges=f['lon_edges'];lat_edges=f['lat_edges']

    maps=[f[f'{k}_field'].reshape(len(lat_edges)-1,len(lon_edges)-1)*100 for k in ('total','genesis','propagation')]

    vmax=max(np.max(np.abs(x)) for x in maps)

    fig=plt.figure(figsize=(12,7.6),constrained_layout=True);gs=fig.add_gridspec(3,2,height_ratios=[1.0,1.0,0.055])

    axes=np.array([[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1])],[fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])]])

    cax=fig.add_subplot(gs[2,0])

    for ax,arr,name in zip(axes.flat[:3],maps,['Total change','Genesis distribution','Conditional track field']):

        ax.text(0.055,0.982,name,transform=ax.transAxes,va='top',ha='left',fontsize=10.5,zorder=10,bbox=dict(facecolor='white',edgecolor='none',alpha=0.75,pad=1.2))

        m=ax.pcolormesh(lon_edges,lat_edges,arr,cmap='RdBu_r',vmin=-vmax,vmax=vmax,shading='flat')

        add_coast(ax,(100,0,180,40));format_geo_axes(ax);ax.set_xticks([100,120,140,160,180]);ax.set_yticks([0,10,20,30,40])

    open_panel_axis(axes[0,1],right=True);axes[0,1].spines['top'].set_visible(True);axes[0,1].spines['left'].set_visible(True)

    axes[0,1].tick_params(axis='y',left=False,labelleft=False,right=True,labelright=False);axes[0,1].yaxis.tick_right();axes[0,1].set_ylabel('')

    axes[1,0].tick_params(axis='x',labelbottom=False)

    fig.colorbar(m,cax=cax,orientation='horizontal',label='Contribution to the 1996–2025 minus 1966–1995 share (percentage points per cell)')

    ax=axes[1,1];x=np.arange(len(s));width=0.34

    g=s.genesis_projection_fraction*100;p=s.propagation_projection_fraction*100

    ax.bar(x-width/2,g,width,color='#4C78A8',label='Genesis distribution')

    ax.bar(x+width/2,p,width,color='#E07B39',label='Conditional track field')

    lo=p-s.propagation_projection_fraction_boot_lo*100;hi=s.propagation_projection_fraction_boot_hi*100-p

    ax.errorbar(x+width/2,p,yerr=[lo,hi],fmt='none',ecolor='0.25',capsize=3,lw=0.9)

    ax.set_xticks(x,[DISPLAY[n] for n in s.index]);ax.set_ylabel('Projection contribution (%)');ax.set_ylim(0,105)

    ax.legend(frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(0.5,-0.16))

    axes[1,1].set_box_aspect(40/80/np.cos(np.deg2rad(20)))

    for panel in axes.flat[:3]:panel.set_xlim(100,180);panel.set_ylim(0,40)

    open_panel_axis(ax,right=True);fig.canvas.draw()

    pos=cax.get_position();cax.set_position([pos.x0,pos.y0-0.004,pos.width,pos.height])

    panel_letters39(axes);save(fig,'Fig04_genesis_conditional_track_field_decomposition',scale=1.10,dpi=320,number=4,title='Genesis-distribution and conditional-track-field decomposition')



# ------------------------------------------------------------------ Figure 6

def figure6():

    direct=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/landfall_latitude_summary.csv')

    direct=direct.loc[direct.start.eq(1966)&direct.end.eq(2025)&direct.metric.eq('mean_lat')].rename(

        columns={'period_difference':'eff','period_ci_low':'lo','period_ci_high':'hi'})

    hd=H('landfall.csv').query("agency=='HKO'").rename(columns={'period_difference':'eff','period_ci_low':'lo','period_ci_high':'hi'})

    direct=pd.concat([direct[['agency','definition','eff','lo','hi']],hd[['agency','definition','eff','lo','hi']]])

    common=pd.read_csv(output_path('results/agency_comparison/common_storm_landfall_four_agency.csv'))  # 1,465 storms common to all four records

    cut=pd.read_csv(ZA/'02_cutpoint_sensitivity/cutpoint_sensitivity.csv')

    cut=cut.loc[cut.scheme.eq('full_record')&cut.metric.eq('first_landfall_latitude_difference'),['cutpoint','agency','effect','p_value']]

    cut=pd.concat([cut,pd.read_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/cutpoint.csv').query("agency=='HKO'")])

    coast=pd.read_csv(ZA/'06_landfall_grouping/coast_grouping_sensitivity.csv')

    coast=coast.loc[coast.threshold_km.eq(50.0)&coast.taiwan_rule.eq('taiwan_north')&coast.denominator.eq('named_only')&coast.event_rule.isin(['all_events','first_any','first_named']),['agency','event_rule','change_percentage_points']]

    coast=pd.concat([coast,H('coast.csv').query("agency=='HKO'")[['agency','event_rule','change_percentage_points']]])

    fig,axes=plt.subplots(2,2,figsize=(11,7.4),constrained_layout=True)

    offsets={'all_events':-0.13,'first_landfall':0.13};top_legends=[]

    for ax,table,e,l,h in [(axes[0,0],direct,'eff','lo','hi'),(axes[0,1],common,'late_minus_early','ci_low','ci_high')]:

        for definition,marker,label in [('all_events','o','All events'),('first_landfall','s','First landfall')]:

            part=table.loc[table.definition.eq(definition)].set_index('agency').loc[AG4]

            x=np.arange(4)+offsets[definition];eff=part[e].to_numpy();low=part[l].to_numpy();high=part[h].to_numpy()

            ax.errorbar(x,eff,yerr=[eff-low,high-eff],fmt=marker,ms=5,capsize=3,color='0.25')

        ax.axhline(0,color='0.5',lw=0.8)

        ax.set(xticks=np.arange(4),xticklabels=['USA','JMA','CMA','HKO'])

        # Marker-only legend entries (the errorbar glyphs were hard to read).

        top_legends.append(ax.legend(handles=[Line2D([0],[0],marker='o',ls='none',color='0.25',ms=5,label='All events'),Line2D([0],[0],marker='s',ls='none',color='0.25',ms=5,label='First landfall')],frameon=False,loc='lower right'))

    axes[0,0].set_ylabel('Latitude change (°)');axes[0,1].set_ylabel('Common-storm latitude change (°)')

    rules=['all_events','first_any','first_named'];rule_labels=['All events','First landfall\n(any coast)','First landfall\n(named coast)']

    for agency in AG4:

        part=coast.loc[coast.agency.eq(agency)].set_index('event_rule').loc[rules]

        axes[1,0].plot(np.arange(3),part.change_percentage_points,marker='o',color=COLORS[agency],label=DISPLAY[agency])

    axes[1,0].axhline(0,color='0.5',lw=0.8)

    axes[1,0].set(xticks=np.arange(3),xticklabels=rule_labels,ylabel='North-share change (percentage points)')

    axes[1,0].legend(frameon=False,ncol=4);axes[1,0].set_ylim(2.0,11.3)  # bottom in the reference layout; top leaves room for the panel letter

    for agency in AG4:

        part=cut.loc[cut.agency.eq(agency)].sort_values('cutpoint')

        axes[1,1].plot(part.cutpoint,part.effect,color=COLORS[agency],lw=1.1,label=DISPLAY[agency])

        sig=part.p_value.lt(0.05)

        axes[1,1].scatter(part.cutpoint[sig],part.effect[sig],s=26,facecolors=COLORS[agency],edgecolors=COLORS[agency],linewidths=0.9,zorder=3)

        axes[1,1].scatter(part.cutpoint[~sig],part.effect[~sig],s=26,facecolors='white',edgecolors=COLORS[agency],linewidths=1.0,zorder=3)

    axes[1,1].axvline(1996,color='0.4',ls='--',lw=0.9);axes[1,1].axhline(0,color='0.5',lw=0.8)

    axes[1,1].set(xlabel='First year of late period',ylabel='First-landfall latitude change (°)');axes[1,1].set_ylim(0.4,1.58)  # bottom in the reference layout; top leaves room for the panel letter

    status=[Line2D([0],[0],marker='o',linestyle='none',markerfacecolor='0.3',markeredgecolor='0.3',markersize=5,label='p < 0.05'),

            Line2D([0],[0],marker='o',linestyle='none',markerfacecolor='white',markeredgecolor='0.3',markersize=5,label='p ≥ 0.05')]

    handles=[Line2D([0],[0],color=COLORS[a],lw=1.2,label=DISPLAY[a]) for a in AG4]

    axes[1,1].legend(handles=handles+status,frameon=False,ncol=3)

    open_panel_axes(axes)

    fig.canvas.draw()

    for ax,legend in zip(axes[0,:],top_legends):

        legend.set_bbox_to_anchor((1.0,12.7/(ax.get_window_extent().height*72.0/fig.dpi)),transform=ax.transAxes)

    panel_letters39(axes)

    save(fig,'Fig06_landfall_latitude_and_grouping',scale=1.27,dpi=600,number=6,title='Landfall latitude and named-coast sensitivity')



# ------------------------------------------------------------------ Figure 7

def figure7():

    """Diagnostic content with consistent figure typography (Arial, bold panel letters)."""

    a=np.load(_LAYOUT_ROOT / 'results/tracks/absolute_track_point_change.npz');change=a['change'];lim=float(np.abs(change).max())

    fig,axs=plt.subplots(1,2,figsize=(12,4.7),gridspec_kw={'width_ratios':[1,1.1]},constrained_layout=True)

    m=axs[0].pcolormesh(a['lon_edges'],a['lat_edges'],change,cmap='RdBu_r',norm=TwoSlopeNorm(0,-lim,lim),shading='flat')

    add_coast(axs[0],(100,0,180,40));format_geo_axes(axs[0]);axs[0].set_xticks([100,120,140,160,180])

    fig.colorbar(m,ax=axs[0],label='1996–2025 minus 1966–1995\n(TS track points per year)',shrink=.8)

    ax=axs[1];ax.set_axis_off()

    boxes=[(.70,'Observed changes','Relative occupancy and absolute passage by latitude band\nTwo-period shoreline-crossing latitude contrast','#dbe8f2'),

           (.38,'Statistical associations','Steering, adjusted-height and PDO covariance\nMonthly sampling decomposition','#e9ecef'),

           (.06,'Processes requiring further tests','Storm-matched steering regimes, vortex structure,\nfine-scale starting positions and storm lifetime','#f2ead8')]

    for y,title,body,color in boxes:

        ax.add_patch(FancyBboxPatch((.03,y),.94,.24,boxstyle='round,pad=.015',facecolor=color,edgecolor='#65727d',linewidth=.7,transform=ax.transAxes))

        ax.text(.5,y+.165,title,ha='center',va='center',fontweight='bold',transform=ax.transAxes,fontsize=10)

        ax.text(.5,y+.07,body,ha='center',va='center',transform=ax.transAxes,fontsize=9,linespacing=1.4)

    panel_letters39(axs)

    save(fig,'Fig07_outlook_schematic',scale=1.10,dpi=320,number=7,title='Absolute track change and levels of evidence')



def figure5():

    """scripts/39 fig5 with two readability additions: a 0.5 m s-1 reference vector and the

    names 'Raw' and 'Detrended' next to the panel letters. Data: archived circulation outputs."""

    fields=np.load(ZV/'circulation_metric_validation/results/wnp_tc_circulation_fields.npz')

    coefs=pd.read_csv(ZV/'circulation_metric_validation/results/wnp_tc_regression_models.csv')

    lon=fields['longitude'];lat=fields['latitude']

    vmax=max(np.max(np.abs(fields['raw_z_eddy_beta'])),np.max(np.abs(fields['detrended_z_eddy_beta'])))

    fig=plt.figure(figsize=(9.15,7.1))

    gs=fig.add_gridspec(2,3,left=0.105,right=0.985,bottom=0.085,top=0.975,width_ratios=[0.045,1.0,1.0],height_ratios=[1.08,0.92],hspace=0.10,wspace=0.10)

    axes=[fig.add_subplot(gs[0,1]),fig.add_subplot(gs[0,2]),fig.add_subplot(gs[1,:])];cax=fig.add_subplot(gs[0,0])

    for ax,scale,name in zip(axes[:2],['raw','detrended'],['Raw','Detrended']):

        m=ax.pcolormesh(lon,lat,fields[f'{scale}_z_eddy_beta'],cmap='RdBu_r',vmin=-vmax,vmax=vmax,shading='auto')

        step=2

        q=ax.quiver(lon[::step],lat[::step],fields[f'{scale}_u_beta'][::step,::step],fields[f'{scale}_v_beta'][::step,::step],color='0.15',scale=18,width=0.0022)

        mask=fields[f'{scale}_u_q']<0.05;yy,xx=np.where(mask)

        if len(xx):ax.scatter(lon[xx],lat[yy],s=3,color='black',alpha=0.45)

        add_coast(ax,(80,0,180,65));format_geo_axes(ax)

        ax.text(0.075,0.985,name,transform=ax.transAxes,va='top',ha='left',fontsize=10,zorder=11,bbox=dict(facecolor='white',edgecolor='none',alpha=0.8,pad=1.2))

        ax.add_patch(plt.Rectangle((0.64,0.012),0.35,0.085,transform=ax.transAxes,facecolor='white',edgecolor='none',alpha=0.85,zorder=8))

        ax.quiverkey(q,0.71,0.055,0.5,r'0.5 m s$^{-1}$',labelpos='E',coordinates='axes',fontproperties={'size':9},zorder=9)

    axes[0].set_xticks([80,100,120,140,160]);axes[1].set_xticks([100,120,140,160,180]);axes[1].tick_params(axis='y',labelleft=False)

    cb=fig.colorbar(m,cax=cax,orientation='vertical');cb.set_label('Regional-background-\nadjusted Z500\nregression\n(m per index SD)',labelpad=7)

    cax.yaxis.set_ticks_position('left');cax.yaxis.set_label_position('left')

    fig.canvas.draw();pos=cax.get_position();pa=axes[0].get_position();cax.set_position([pos.x0-0.045,pa.y0,pos.width,pa.height])

    keep=coefs[(coefs.term!='const')&coefs.model.isin(['eddy_wnpsh_oni_pdo','steering_oni_pdo'])]

    keep=keep[keep.term.isin(['eddy_wnpsh_mean_m','corridor_u_steer_ms','jas_oni','pdo'])]

    order=[('raw','eddy_wnpsh_mean_m','Regional-background-adjusted\nZ500 index'),('detrended','eddy_wnpsh_mean_m','Regional-background-adjusted\nZ500 index, detrended'),

           ('raw','corridor_u_steer_ms','Zonal steering flow'),('detrended','corridor_u_steer_ms','Zonal steering flow, detrended'),

           ('raw','pdo','PDO'),('detrended','pdo','PDO, detrended')]

    rows=[]

    for scale,term,label in order:

        r=keep[(keep.timescale==scale)&(keep.term==term)].iloc[0];rows.append((label,r.coefficient_standardized,r.ci_low,r.ci_high,r.q_bh_within_timescale))

    y=np.arange(len(rows));beta=np.array([r[1] for r in rows]);lo=beta-np.array([r[2] for r in rows]);hi=np.array([r[3] for r in rows])-beta

    axes[2].errorbar(beta,y,xerr=[lo,hi],fmt='o',color='#365E8D',ecolor='0.25',capsize=3)

    axes[2].axvline(0,color='0.45',lw=0.8);axes[2].set_yticks(y,[r[0] for r in rows]);axes[2].invert_yaxis();axes[2].set_xlabel('Standardized coefficient (HAC 95% CI)')

    for yi,r in enumerate(rows):axes[2].text(r[3]+0.004,yi,f'q={r[4]:.3f}',va='center',fontsize=10)

    open_panel_axis(axes[2]);panel_letters39(axes)

    save(fig,'Fig05_circulation_context',scale=1.10,dpi=320,number=5,title='Circulation context of the annual redistribution index')



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    for k in (sys.argv[1:] or ['2','3','4','5','6']):globals()[f'figure{k}']()
