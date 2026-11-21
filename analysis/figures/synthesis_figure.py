"""Render Figure 6, the pathway-change schematic."""
from __future__ import annotations

from analysis.figures.figure_raster import render_png, resize_terrain, apply_tnr
import argparse
from pathlib import Path

import numpy as np

from core.paths import PROJECT, output_path

EXTENT=(105.0,175.0,0.0,50.0)

BASE_MAP_W=429.0*80/50*np.cos(np.deg2rad(25))

MAP_W=BASE_MAP_W*70/80*1.18

MAP_H=round(BASE_MAP_W*50/80/np.cos(np.deg2rad(25))*1.18,1)

MAP_TOP=6.0+(243.1-MAP_H)/2

WIDTH=3600

ARROW=(0.29,0.29,0.29)

LANDFALL=(0.184,0.361,0.541)

FRAME=(0.729,0.765,0.788)

# Colours used by the author-approved cropped-map source.
WEST,WEST_TEXT=(0.706,0.396,0.247),'#8c4a2b'
RECURVE,RECURVE_TEXT=(0.184,0.431,0.561),'#22536b'

def draw(out):
    output = Path(out).expanduser()
    if not output.is_absolute():
        output = PROJECT / output
    output = output.resolve()
    if not output.is_relative_to(PROJECT.resolve()):
        raise ValueError('Output must remain within the project.')
    if output.suffix:
        raise ValueError('Pass an output stem without a file extension.')
    if any(Path(str(output) + suffix).exists() for suffix in ('.png', '.pdf')):
        raise FileExistsError('Refusing to overwrite an existing figure.')
    output.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    import matplotlib.patheffects as pe

    from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Rectangle

    from matplotlib.colors import TwoSlopeNorm

    import geopandas as gpd

    import analysis.figures.main_figures as mf   # font set-up and geographic tick format

    FONT=mf.FONT

    W,H=672.0+(BASE_MAP_W-324.3)+(667.6-364.7)*0.30+BASE_MAP_W*70/80*0.18,630.0                     # space below the map for its longitude labels and colour bar

    fig=plt.figure(figsize=(W/72,H/72))

    fx=lambda x:x/W;fy=lambda y:1-y/H  # PDF points (origin top-left) to figure fraction



    # ---- (a) map frame: x 26.3-350.6 pt, y 15.75-239.35 pt; colour bar below the map
    ax=fig.add_axes([fx(26.3),fy(4.0+MAP_H),MAP_W/W,MAP_H/H])
    # Same ETOPO1 colours, hillshade and blending used for the Figure 5 maps.
    import xarray as xr
    from matplotlib.colors import LightSource, LinearSegmentedColormap
    with xr.open_dataset(PROJECT / 'data/derived/etopo1_100-180E_0-50N_2min.nc') as dataset:
        terrain=dataset.altitude.sel(latitude=slice(EXTENT[2],EXTENT[3]),longitude=slice(EXTENT[0],EXTENT[1]))
        terrain_extent=[float(terrain.longitude.min()),float(terrain.longitude.max()),float(terrain.latitude.min()),float(terrain.latitude.max())]
        z=terrain.values.astype(float)
        if terrain.latitude.values[0] < terrain.latitude.values[-1]: z=z[::-1]
    light=LightSource(azdeg=315,altdeg=40)
    land_cmap=LinearSegmentedColormap.from_list('land',['#C9D3AE','#D9D2A6','#C7B48C','#A78F72','#8C7A68','#F1EEEA'])
    sea_cmap=LinearSegmentedColormap.from_list('sea',['#9FB7CB','#BFD0DE','#DCE7EF','#EEF3F7'])
    land_rgb=light.shade(np.clip(z,0,4500),cmap=land_cmap,vmin=0,vmax=4500,vert_exag=0.04,blend_mode='soft')
    sea_rgb=light.shade(np.clip(z,-7000,0),cmap=sea_cmap,vmin=-7000,vmax=0,vert_exag=0.006,blend_mode='soft')
    relief=np.where((z>0)[...,None],0.45*land_rgb+0.55,0.6*sea_rgb+0.4)
    ax.imshow(relief,extent=terrain_extent,origin='upper',aspect='auto',interpolation='bilinear',zorder=0)
    a=np.load(PROJECT / 'results/tracks/absolute_track_point_change.npz')
    change=np.ma.masked_where(np.abs(a['change'])<0.25,a['change'])
    lim=2.0;ALPHA=0.6
    m=ax.pcolormesh(a['lon_edges'],a['lat_edges'],change,cmap='RdBu_r',norm=TwoSlopeNorm(0,-lim,lim),shading='flat',alpha=ALPHA,zorder=1,linewidth=0,rasterized=True)
    land=gpd.read_file(PROJECT/'data/raw/GSHHG/GSHHS_shp/i/GSHHS_i_L1.shp',bbox=(EXTENT[0],EXTENT[2],EXTENT[1],EXTENT[3]))
    land.boundary.plot(ax=ax,color='0.40',linewidth=0.3,zorder=2)
    ax.set_xlim(*EXTENT[:2]);ax.set_ylim(*EXTENT[2:])
    ax.set_xticks([110,130,150,170]);ax.set_yticks([0,10,20,30,40,50])
    mf.format_geo_axes(ax);ax.tick_params(labelsize=13.3,width=0.8,length=2.4)
    ax.tick_params(axis='x', top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.spines['top'].set_visible(True)
    ax.spines['top'].set_color('0.15')
    ax.spines['top'].set_zorder(20)
    for s in ax.spines.values():s.set_linewidth(0.7)
    lon=lambda x:100+(x-26.3)/4.054;lat=lambda y:(249.1-y)/4.862   # reference PDF points to degrees
    halo=[pe.withStroke(linewidth=2.6,foreground='white',alpha=0.85)]
    ax.add_patch(FancyArrowPatch((lon(280.0),lat(185.2)),(lon(149.3),lat(168.2)),connectionstyle='arc3,rad=0.12',
        arrowstyle='-|>,head_length=4.2,head_width=2.6',color=WEST,lw=3.4,zorder=4,shrinkA=0,shrinkB=0))
    ax.add_patch(FancyArrowPatch((lon(264.0),lat(155.8)),(lon(205.3),lat(66.0)),connectionstyle='arc3,rad=-0.28',
        arrowstyle='-|>,head_length=4.0,head_width=2.4',color=RECURVE,lw=4.4,zorder=4,shrinkA=0,shrinkB=0))
    ax.add_patch(FancyArrowPatch((lon(129.6),lat(128.0)),(lon(129.6),lat(90.0)),arrowstyle='-|>,head_length=4.2,head_width=2.6',
        color=LANDFALL,lw=3.4,zorder=4,shrinkA=0,shrinkB=0))
    t=dict(fontsize=12.8,zorder=5,path_effects=halo)
    ax.text(lon(237.1),lat(212.1),'Low-latitude westward corridor',color=WEST_TEXT,ha='center',va='center',**t)
    ax.text(lon(237.1),lat(227.5),'fewer passages',color=WEST_TEXT,ha='center',va='center',style='italic',**t)
    ax.text(lon(263.0),lat(95.4),'Recurving corridor',color=RECURVE_TEXT,ha='center',va='center',**t)
    ax.text(lon(263.0),lat(110.9),'larger relative share',color=RECURVE_TEXT,ha='center',va='center',style='italic',**t)
    ax.text(lon(50.0),lat(107.5),'Poleward shift',color=LANDFALL,ha='left',va='center',**t)
    ax.text(lon(50.0),lat(123.0),'in landfall',color=LANDFALL,ha='left',va='center',**t)
    ax.text(lon(50.0),lat(138.5),'latitude',color=LANDFALL,ha='left',va='center',**t)
    cbar_width=0.76*MAP_W
    cax=fig.add_axes([fx(26.3+(MAP_W-cbar_width)/2),fy(4.0+MAP_H+40.0),cbar_width/W,14.25/H])
    cb=fig.colorbar(m,cax=cax,orientation='horizontal',extend='both',ticks=[-2,-1,0,1,2])
    cb.solids.set_alpha(ALPHA);cb.outline.set_linewidth(0.5)
    cax.tick_params(labelsize=13.3,length=1.8,width=0.5,pad=1.2)
    cb.ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,_:f'{v:g}'.replace('-','−')))
    cax.set_xlabel('Change in annual track points per cell',fontsize=14.0,labelpad=2.5,linespacing=1.2)
    ax.text(0.012,0.985,'a',transform=ax.transAxes,ha='left',va='top',fontsize=16.7,fontweight='bold',zorder=8)

    # ---- (b) shared result, statistical analyses, interpretation, untested processes, and expressions
    BS,TS=12.5,13.2
    SB,ST=12.1,12.7                      # body and title sizes in the three narrow boxes
    YB=4.0+MAP_H+114.0  # extend panel b to the bottom of the figure
    bx=fig.add_axes([fx(364.7+(MAP_W-324.3)),fy(YB),(1.19*(1.30*(667.6-364.7)+(BASE_MAP_W-BASE_MAP_W*70/80)))/W,(YB-4.0)/H])
    bx.set_xlim(364.7,667.6);bx.set_ylim(455.0,4.0);bx.set_axis_off()
    def box(x0,y0,x1,y1,kind,items):
        styles={
            'observed':('#315b83','#234360','white','solid',1.0),
            'partition':('#eae5d6','#97875f','#39362b','solid',1.0),
            'interpretation':('#f4f7fa','#315b83','#234360',(0,(4,2.5)),1.2),
            'untested':('#ffffff','#7d8285','#4a4f52',(0,(1,2)),1.3),
            'expression':('#dce6ef','#7795ad','#294354','solid',1.0),
        }
        fc,ec,tc,ls,lw=styles[kind]
        bx.add_patch(FancyBboxPatch((x0,y0),x1-x0,y1-y0,
            boxstyle='round,pad=0,rounding_size=4.5',facecolor=fc,edgecolor=ec,lw=lw,linestyle=ls))
        for index,(text,offset,size,weight) in enumerate(items):
            # Redistribute text vertically to fit the enlarged font.
            if len(items) > 1:
                offset = 12 + index * (y1-y0-24)/(len(items)-1)
            bx.text((x0+x1)/2,y0+offset,text,ha='center',va='center',
                fontsize=size,fontweight=weight,color=tc)
    def segment(x0,y0,x1,y1,ls='solid'):
        bx.plot([x0,x1],[y0,y1],color='#4b5357',lw=1.0,linestyle=ls,solid_capstyle='butt',zorder=1)
    def arrow(x0,y0,x1,y1,style='-|>',ls='solid'):
        bx.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle=style,
            mutation_scale=10,lw=1.0,color='#4b5357',linestyle=ls,shrinkA=0,shrinkB=0,zorder=2))
    X0,X1=368.7,663.6
    GAP=6.0;CW=(X1-X0-2*GAP)/3
    c1=(X0,X0+CW);c2=(X0+CW+GAP,X0+2*CW+GAP);c3=(X1-CW,X1)
    m1=sum(c1)/2;m2=sum(c2)/2;m3=sum(c3)/2;m12=(c1[0]+c2[1])/2
    L,R=436.45,587.85                  # centres of the two expression boxes
    # row 1: shared observed result
    box(X0,12,X1,57,'observed',[
        ('Shared shift toward higher-latitude pathways',16,12.8,'bold'),
        ('Relative track composition across four agencies',32,BS,'normal')])
    # row 2: two partitions (arrows) and one association (double-headed line, no causal direction)
    Y2,Y2b=88,225
    segment(m12,57,m12,66);segment(m1,66,m2,66)
    arrow(m1,66,m1,Y2-2);arrow(m2,66,m2,Y2-2)
    bx.text(m12,78,'Partitions',ha='center',va='center',fontsize=BS,color='#555d61')
    arrow(m3,57,m3,Y2-2,style='<|-|>')
    bx.text(m3-6,74,'Association',ha='right',va='center',fontsize=BS,color='#555d61')
    box(c1[0],Y2,c1[1],Y2b,'partition',[
        ('Genesis-bin',13,ST,'bold'),('partition',26,ST,'bold'),
        ('74%–80%',50,14,'bold'),
        ('of projected',73,SB,'normal'),('contrast in',87,SB,'normal'),
        ('conditional',101,SB,'normal'),('track fields',114,SB,'normal')])
    box(c2[0],Y2,c2[1],Y2b,'partition',[
        ('Six-hourly',13,ST,'bold'),('steering',26,ST,'bold'),
        ('Primary record,',45,SB,'normal'),('June–October',59,SB,'normal'),
        ('Path changes lie',78,SB,'normal'),('mainly within',92,SB,'normal'),
        ('defined regimes',106,SB,'normal')])
    box(c3[0],Y2,c3[1],Y2b,'partition',[
        ('Circulation',13,ST,'bold'),('association',26,ST,'bold'),
        ('PDO and zonal',50,SB,'normal'),('steering covary',64,SB,'normal'),
        ('with the index;',78,SB,'normal'),('causal order',97,SB,'normal'),
        ('not tested',111,SB,'normal')])
    # row 3: spatial interpretation of the partitions; untested processes below the association
    Y3,Y3b=246,321
    segment(m1,Y2b,m1,Y2b+2);segment(m2,Y2b,m2,Y2b+2);segment(m1,Y2b+2,m2,Y2b+2)
    arrow(m12,Y2b+2,m12,Y3-2)
    bx.text(m12+6,Y2b+11,'Interpretation',ha='left',va='center',fontsize=BS,color='#555d61')
    box(c1[0],Y3,c2[1],Y3b,'interpretation',[
        ('Different routes can share',14,TS,'bold'),
        ('broad genesis regions',28,TS,'bold'),
        ('Regional track shares can shift',44,BS,'normal'),
        ('in opposite directions',57,BS,'normal')])
    arrow(m3,Y2b,m3,Y3-2,ls=(0,(2,2)))
    box(c3[0],Y3,c3[1],Y3b,'untested',[
        ('Untested',14,ST,'bold'),('processes',27,ST,'bold'),
        ('Circulation-to-',44,SB,'normal'),('track links',57,SB,'normal')])
    # row 4: along-track expressions of the shared shift
    Y4,Y4b=355,451
    segment(m12,Y3b,m12,Y3b+13);segment(L,Y3b+13,R,Y3b+13)
    arrow(L,Y3b+13,L,Y4-2);arrow(R,Y3b+13,R,Y4-2)
    bx.text((L+R)/2,Y3b+23,'Along-track expressions',ha='center',va='center',fontsize=BS,color='#555d61')
    box(368.7,Y4,504.2,Y4b,'expression',[
        ('Coastal',15,TS,'bold'),
        ('displacement',30,TS,'bold'),
        ('Higher mean',32,BS,'normal'),
        ('first-landfall latitude',48,BS,'normal'),
        ('in the later period',63,BS,'normal')])
    box(512.1,Y4,663.6,Y4b,'expression',[
        ('Upper-ocean',12,TS,'bold'),
        ('exposure',24,TS,'bold'),
        ('Path and ocean changes',40,BS,'normal'),
        ('act in opposite directions;',56,BS,'normal'),
        ('net change unresolved',72,BS,'normal')])
    bx.text(355.8,6,'b',ha='left',va='top',fontsize=16.7,fontweight='bold')
    for t in fig.findobj(match=matplotlib.text.Text):
        t.set_fontfamily(FONT)
        if t.get_visible() and t.get_fontsize()>0:
            t.set_fontsize(t.get_fontsize()*1.90)

    meta={'Title':'Pathway change, statistical analyses, and untested processes','Subject':'Journal of Geophysical Research: Atmospheres Figure 6','Creator':'Matplotlib'}

    stem=str(output)

    apply_tnr(fig)
    fig.savefig(output_path(stem+'.pdf'),bbox_inches='tight',pad_inches=0.02,metadata=meta,dpi=600)
    pdf = output.with_suffix('.pdf')
    resize_terrain(pdf, (4777, 4218), 0.34, 1)
    render_png(pdf)

    plt.close(fig)
    print('Figure saved:',stem)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true', help='Render the figure.')
    parser.add_argument('--output', required=True, help='New project-local output stem without an extension.')
    args = parser.parse_args()
    if not args.execute:
        parser.error('Add --execute to render the figure.')
    draw(args.output)


if __name__ == '__main__':
    main()
