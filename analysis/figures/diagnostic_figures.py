'Diagnostic figures for track counts, circulation and ocean exposure.'
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

from core.paths import PROJECT,output_path,write_json

import argparse

import numpy as np

import pandas as pd

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib.patches import FancyBboxPatch,Rectangle

from matplotlib.lines import Line2D

from matplotlib.colors import TwoSlopeNorm

import geopandas as gpd





from analysis.figures.shared_figure_style import COLORS as COLOURS, LABELS as LABEL, PERIODS, geo_ticks, FONT

EARLY,LATE=PERIODS['early'][0],PERIODS['late'][0]

_coast=None



# Canvas dimensions in pixels for the supplementary figures.

CANVASES={

    'FigS14_absolute_passage_recurvature':(3184,1018),

    'FigS15_bam_baseline_skill':(3632,1294),

    'FigS16_steering_layer_regressions':(3332,2269),

    'FigS17_circulation_geometry_diagnostics':(3379,1153),

    'FigS18_hko_common_window':(3029,1541),

    'FigS19_d26_parallel_exposure':(3334,1204),

}

plt.rcParams.update({'font.size':11.5,'axes.titlesize':11.5,'axes.labelsize':11.5,

                     'xtick.labelsize':10.5,'ytick.labelsize':10.5,'legend.fontsize':10.5})



def canvas(name):

    # Guard against floating-point truncation of an exact integer pixel count.

    return tuple((v+.01)/300 for v in CANVASES[name])



def panel(ax,letter,title=''):

    ax.text(.015,.98,letter,transform=ax.transAxes,ha='left',va='top',

            fontsize=12.5,fontweight='bold',zorder=12,

            bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=.6))

    if title:ax.set_title(title,pad=7)

def coast(ax,extent=(100,180,0,40)):

    global _coast

    if _coast is None:_coast=gpd.read_file(PROJECT/'data/raw/GSHHG/GSHHS_shp/i/GSHHS_i_L1.shp',bbox=(80,0,180,65))

    _coast.plot(ax=ax,facecolor='#E6E4DF',edgecolor='#505050',linewidth=.35,zorder=3)

    ax.set_xlim(extent[0],extent[1]);ax.set_ylim(extent[2],extent[3])

    # Map scale true at the central latitude of the panel (same convention as Figure 1).

    ax.set_aspect(1/np.cos(np.deg2rad((extent[2]+extent[3])/2)),adjustable='box')

    # Degree-labelled ticks as in the main figures and the other SI maps.

    geo_ticks(ax,list(range(int(extent[0]),int(extent[1])+1,20)),list(range(int(extent[2]),int(extent[3])+1,10)))



def save(fig,name):

    # All typography is set before layout. No post-layout font enlargement.

    fig.savefig(output_path('figures/'+name+'.png'),dpi=300)

    fig.savefig(output_path('figures/'+name+'.pdf'),metadata={

        'Creator':'Matplotlib','Author':'Daizhao Huang; Yuanfu Yue'})

    plt.close(fig);print('Figure',name,flush=True)



def absolute():

    a=pd.read_csv(_LAYOUT_ROOT / 'results/tracks/absolute_track_periods.csv');r=pd.read_csv(_LAYOUT_ROOT / 'results/tracks/recurvature_summary.csv')

    agencies=['PRIMARY','USA','TOKYO','CMA','HKO'];x=np.arange(5);fig,axs=plt.subplots(1,3,figsize=canvas('FigS14_absolute_passage_recurvature'),layout='constrained')

    for ax,band,letter,title in zip(axs[:2],['0-20','20-40'],'ab',['0°–20°N','20°–40°N']):

        for k,(p,c) in enumerate([('early',EARLY),('late',LATE)]):

            v=a.loc[a.latitude_band.eq(band)&a.period.eq(p)].set_index('agency').loc[agencies,'ts_track_points']

            ax.bar(x+(k-.5)*.34,v,width=.32,color=c,label='1966–1995' if k==0 else '1996–2025')

        ax.set_xticks(x,[LABEL[k] for k in agencies]);ax.set_ylabel('Tropical-storm track points per year');panel(ax,letter,title)

    for k,(p,c) in enumerate([('early',EARLY),('late',LATE)]):

        v=r.loc[r.period.eq(p)].set_index('agency').loc[agencies,'recurving_per_year']

        axs[2].bar(x+(k-.5)*.34,v,width=.32,color=c)

    axs[2].set_xticks(x,[LABEL[k] for k in agencies]);axs[2].set_ylabel('Recurving storms per year');panel(axs[2],'c','Recurving storms')

    handles,labels=axs[0].get_legend_handles_labels()

    fig.legend(handles,labels,loc='outside lower center',ncol=2,frameon=False)

    save(fig,'FigS14_absolute_passage_recurvature')



def bam():

    # Primary calibration (elapsed time between successive reports), as in Table S9.

    a=np.load(_LAYOUT_ROOT / 'results/models/bam_sensitivity_annual.npz')

    o=a['observed'];obs=[o[:30].mean(0).reshape(16,32),o[30:].mean(0).reshape(16,32)]

    sim=[a['850-300_actual_dt_max12h_'+k].mean(0).reshape(16,32) for k in ['early','late']]

    fig,axs=plt.subplots(2,3,figsize=canvas('FigS15_bam_baseline_skill'))

    fig.subplots_adjust(left=.045,right=.985,bottom=.22,top=.90,wspace=.18,hspace=.24)

    for row,p in enumerate(['1966–1995','1996–2025']):

        for col,(field,title) in enumerate([(obs[row],'Observed'),(sim[row],'Beta-advection model'),(sim[row]-obs[row],'Model minus observation')]):

            if col<2:im=axs[row,col].imshow(field*100,origin='lower',extent=[100,180,0,40],cmap='Blues',vmin=0,vmax=max(o.max() for o in obs+sim)*100,aspect='auto')

            else:

                lim=max(abs(sim[k]-obs[k]).max() for k in [0,1])*100

                im=axs[row,col].imshow(field*100,origin='lower',extent=[100,180,0,40],cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')

            coast(axs[row,col]);panel(axs[row,col],chr(97+row*3+col),title if row==0 else '')

            axs[row,col].text(.97,.98,p,transform=axs[row,col].transAxes,

                ha='right',va='top',fontsize=10.5,zorder=12,

                bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1))

            if col==1:dens=im

            if col==2:diff=im

    # One colour bar for the density panels and one for the differences.

    fig.colorbar(dens,cax=fig.add_axes([.075,.105,.51,.025]),orientation='horizontal',label='Relative track composition (%)')

    fig.colorbar(diff,cax=fig.add_axes([.715,.105,.24,.025]),orientation='horizontal',label='Difference (percentage points)')

    save(fig,'FigS15_bam_baseline_skill')



def steering():

    a=np.load(_LAYOUT_ROOT / 'results/circulation/steering_field_regressions.npz');lon=a['longitude'];lat=a['latitude'];X,Y=np.meshgrid(lon,lat)

    fig,axs=plt.subplots(2,2,figsize=canvas('FigS16_steering_layer_regressions'))

    fig.subplots_adjust(left=.06,right=.84,bottom=.065,top=.94,wspace=.17,hspace=.19)

    vmax=max(abs(a[f'{layer}_{scale}_z_eddy_beta']).max() for layer in ['300','500'] for scale in ['raw','detrended'])

    for row,layer in enumerate(['300','500']):

        for col,scale in enumerate(['raw','detrended']):

            ax=axs[row,col]

            im=ax.pcolormesh(lon,lat,a[f'{layer}_{scale}_z_eddy_beta'],cmap='RdBu_r',vmin=-vmax,vmax=vmax,shading='auto')

            q=ax.quiver(X[::2,::2],Y[::2,::2],a[f'{layer}_{scale}_u_beta'][::2,::2],a[f'{layer}_{scale}_v_beta'][::2,::2],scale=12,color='#202020',width=.003)

            ax.add_patch(Rectangle((.73,.015),.26,.08,transform=ax.transAxes,facecolor='white',edgecolor='none',zorder=9))

            ax.quiverkey(q,.80,.055,.5,r'0.5 m s$^{-1}$',labelpos='E',coordinates='axes',zorder=10)

            coast(ax,(80,180,0,65));panel(ax,chr(97+row*2+col),f'850–{layer} hPa, {scale}')

    fig.colorbar(im,cax=fig.add_axes([.88,.15,.018,.70]),label='Regional-background-adjusted\n500-hPa height (m per index SD)')

    save(fig,'FigS16_steering_layer_regressions')

    fig,axs=plt.subplots(1,2,figsize=canvas('FigS17_circulation_geometry_diagnostics'))

    axs[0].set_position([.055,.13,.36,.70]);axs[1].set_position([.625,.13,.36,.70])

    zlim=float(np.nanmax(abs(a['period_change_z_eddy'])))

    im=axs[0].pcolormesh(lon,lat,a['period_change_z_eddy'],cmap='RdBu_r',vmin=-zlim,vmax=zlim,shading='auto')

    for p,c in [('early',EARLY),('late',LATE)]:

        cs=axs[0].contour(lon,lat,a['z500_'+p],levels=[5880],colors=[c],linewidths=1.2,zorder=4)

        axs[0].clabel(cs,fmt={5880:'5880 gpm'},fontsize=10)

        axs[1].contour(lon,lat,a['u500_'+p],levels=[0],colors=[c],linewidths=1.2,zorder=4)

    q=axs[1].quiver(X[::2,::2],Y[::2,::2],a['period_change_u300'][::2,::2],a['period_change_v300'][::2,::2],scale=15,width=.003)

    axs[1].add_patch(Rectangle((.71,.875),.28,.08,transform=axs[1].transAxes,facecolor='white',edgecolor='none',zorder=9))

    axs[1].quiverkey(q,.78,.92,.5,r'0.5 m s$^{-1}$',labelpos='E',coordinates='axes',zorder=10)

    for ax in axs:coast(ax,(80,180,0,65))

    panel(axs[0],'a','Height change and\n5880-gpm contours')

    panel(axs[1],'b','500-hPa zero zonal-wind lines\nand steering change')

    fig.colorbar(im,cax=fig.add_axes([.445,.15,.017,.64]),label='1996–2025 minus 1966–1995\nbackground-adjusted height (m)')

    axs[1].legend(handles=[Line2D([0],[0],color=EARLY,label='1966–1995'),Line2D([0],[0],color=LATE,label='1996–2025')],loc='upper left',bbox_to_anchor=(.07,1),frameon=True,facecolor='white',framealpha=1,edgecolor='white',fontsize=10)

    save(fig,'FigS17_circulation_geometry_diagnostics')



def hko():

    a=np.load(_LAYOUT_ROOT / 'results/tracks/hko_track_fields.npz');fig,axs=plt.subplots(2,2,figsize=canvas('FigS18_hko_common_window'))

    fig.subplots_adjust(left=.06,right=.84,bottom=.085,top=.94,wspace=.20,hspace=.28)

    lim=max(abs(a['native_all_'+k]).max() for k in ['USA','TOKYO','CMA','HKO'])*100

    for ax,k,label in zip(axs.flat,['USA','TOKYO','CMA','HKO'],'abcd'):

        im=ax.imshow(a['native_all_'+k]*100,origin='lower',extent=[100,180,0,40],cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')

        coast(ax);panel(ax,label,LABEL[k])

    fig.colorbar(im,cax=fig.add_axes([.88,.15,.02,.70]),label='Relative track-composition change\n(percentage points)')

    save(fig,'FigS18_hko_common_window')



def evidence():

    from analysis.shared.agency_catalogues import primary_ts, domain

    t=domain(primary_ts());xe=np.arange(100,182.5,2.5);ye=np.arange(0,42.5,2.5)

    counts=np.array([np.histogram2d(t.loc[t.season.eq(y),'lon'],t.loc[t.season.eq(y),'lat'],bins=[xe,ye])[0].T for y in range(1966,2026)])

    change=counts[30:].mean(0)-counts[:30].mean(0)

    np.savez(_LAYOUT_ROOT / 'results/tracks/absolute_track_point_change.npz',change=change,early=counts[:30].mean(0),late=counts[30:].mean(0),lon_edges=xe,lat_edges=ye)







def d26():

    s=pd.read_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_parallel_summary.csv');b=pd.read_csv(_LAYOUT_ROOT / 'results/ocean/d26/d26_parallel_bootstrap.csv')

    terms=['delta_exposure','redistribution_component','ocean_component','covariability_component'];labels=['Net change','Path–season','Ocean field','Annual covariability']

    fig,axs=plt.subplots(1,2,figsize=canvas('FigS19_d26_parallel_exposure'),constrained_layout=True)

    for ax,var,unit,div in zip(axs,['D26','OHC300_same_support'],['m',r'GJ m$^{-2}$'],[1,1e9]):

        row=s.loc[s.variable.eq(var)&s.land_treatment.eq('ocean_points_only')&s.comparison.eq('1966_1995_vs_1996_2025')].iloc[0]

        intervals=b.loc[b.variable.eq(var)&b.land_treatment.eq('ocean_points_only')&b.comparison.eq('1966_1995_vs_1996_2025')].set_index('component').loc[terms]

        values=np.array([row[t] for t in terms])/div;lo=intervals.q025.to_numpy()/div;hi=intervals.q975.to_numpy()/div

        for k,(value,l,h) in enumerate(zip(values,lo,hi)):

            ax.plot([l,h],[k,k],color='#365E8D',lw=1.4);ax.plot(value,k,'o',color='#365E8D')

        ax.axvline(0,color='.5',lw=.7);ax.set_yticks(range(4),labels);ax.invert_yaxis();ax.set_xlabel(f'1996–2025 minus 1966–1995 ({unit})');panel(ax,'ab'[0 if var=='D26' else 1],var.replace('_same_support',''))

    save(fig,'FigS19_d26_parallel_exposure')



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    parser=argparse.ArgumentParser();parser.add_argument('--only',choices=['all','d26','counts_and_circulation','s15','maps','si'],default='all');args=parser.parse_args()

    if args.only=='counts_and_circulation':absolute();steering();evidence()

    if args.only=='all':absolute();bam();steering();hko();evidence()

    if args.only=='s15':bam()

    if args.only=='maps':bam();steering();hko()

    if args.only=='si':absolute();bam();steering();hko();d26()

    if args.only in ['all','d26'] and _LAYOUT_ROOT / 'results/ocean/d26/d26_parallel_summary.csv'.exists():d26()
