'Synthesis figure.'
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import sys,time

import numpy as np

from PIL import Image

import analysis.figures.relief_basemap as mb

from core.paths import PROJECT,output_path,write_json



EXTENT=(100.0,180.0,0.0,50.0)

BASE=_LAYOUT_ROOT / 'data/figure_inputs/mff_relief_z7_fig07.png'

# Panel geometry of the reference figure (PDF points): map frame 324.3 x 243.1 pt.

MAP_W=324.3

MAP_H=round(MAP_W*50/80/np.cos(np.deg2rad(25)),1)   # 223.6 pt: scale true at 25°N, as in Figure 1

MAP_TOP=6.0+(243.1-MAP_H)/2

WIDTH=3600



def tiles():

    x0,y0=mb.tile_xy(EXTENT[0],EXTENT[3]);x1,y1=mb.tile_xy(EXTENT[1]-1e-9,EXTENT[2]+1e-9)

    return [(x,y) for x in range(int(x0),int(x1)+1) for y in range(int(y0),int(y1)+1)]



def download(budget=150):

    from concurrent.futures import ThreadPoolExecutor

    mb.TILE_DIR.mkdir(parents=True,exist_ok=True);todo=tiles();start=time.monotonic();status=[]

    with ThreadPoolExecutor(8) as pool:

        for i in range(0,len(todo),32):

            if time.monotonic()-start>budget:break

            status+=list(pool.map(mb.fetch,todo[i:i+32]))

    have=sum((mb.TILE_DIR/f'{x}_{y}.jpg').exists() for x,y in todo)

    print(f'tiles {have}/{len(todo)}',{s:status.count(s) for s in set(status)})



def build():

    from scipy.ndimage import map_coordinates

    from analysis.figures.study_region_figure import recolour

    todo=tiles();xs=sorted({x for x,_ in todo});ys=sorted({y for _,y in todo})

    mosaic=np.zeros((len(ys)*256,len(xs)*256,3),np.uint8)

    for x,y in todo:

        mosaic[(y-ys[0])*256:(y-ys[0]+1)*256,(x-xs[0])*256:(x-xs[0]+1)*256]=np.array(Image.open(mb.TILE_DIR/f'{x}_{y}.jpg').convert('RGB'))

    height=round(WIDTH*MAP_H/MAP_W)

    lon=EXTENT[0]+(np.arange(WIDTH)+0.5)/WIDTH*(EXTENT[1]-EXTENT[0])

    lat=EXTENT[3]-(np.arange(height)+0.5)/height*(EXTENT[3]-EXTENT[2])

    px=(lon+180)/360*2**mb.ZOOM*256-xs[0]*256

    py=(1-np.arcsinh(np.tan(np.radians(lat)))/np.pi)/2*2**mb.ZOOM*256-ys[0]*256

    PY,PX=np.meshgrid(py,px,indexing='ij')

    rgb=np.stack([map_coordinates(mosaic[...,c].astype(float),[PY-0.5,PX-0.5],order=1,mode='nearest') for c in range(3)],-1).clip(0,255).astype(np.uint8)

    # Figure 1 uses sigma 14 px at 0.012 deg px-1; the same ground scale here is about 8 px.

    grey,sea=recolour(rgb,sigma=8,post=0.5)

    Image.fromarray(grey).save(BASE)

    write_json('results/diagnostics/fig07_basemap.json',dict(status='PASS',source=mb.URL.replace('{z}',str(mb.ZOOM)),zoom=mb.ZOOM,tiles=len(todo),

        output=str(BASE.relative_to(output_path('.'))),size=[WIDTH,height],extent=EXTENT,recolour='study_region_figure.recolour(sigma=8, post=0.5)',sea_fraction=round(float(sea),4)))

    print('basemap',WIDTH,height)



# ------------------------------------------------------------------ drawing

# Colours of the reference figure.

OBS_FILL,OBS_EDGE=(0.184,0.361,0.541),(0.125,0.251,0.373)

ASSOC_FILL,ASSOC_EDGE,ASSOC_TEXT=(0.91,0.886,0.831),(0.541,0.478,0.333),'#3a3320'

OPEN_EDGE,OPEN_TEXT=(0.42,0.42,0.42),'#333333'

ARROW=(0.29,0.29,0.29)

WEST,WEST_TEXT=(0.706,0.396,0.247),'#8c4a2b'

RECURVE,RECURVE_TEXT=(0.184,0.431,0.561),'#22536b'

LANDFALL=(0.184,0.361,0.541)

FRAME=(0.729,0.765,0.788)



def draw(variant='overlay',out=None):

    import matplotlib

    matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    import matplotlib.patheffects as pe

    from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Rectangle

    from matplotlib.colors import TwoSlopeNorm

    import geopandas as gpd

    import analysis.figures.main_figures as mf   # font set-up and geographic tick format

    FONT=mf.FONT

    W,H=672.0,292.0                     # reference figure 672 x 264 pt; 28 pt added for the legend row

    fig=plt.figure(figsize=(W/72,H/72))

    fx=lambda x:x/W;fy=lambda y:1-y/H  # PDF points (origin top-left) to figure fraction



    # ---- (a) map, same frame in the reference layout: x 26.3-350.6 pt, y 6.0-249.1 pt

    ax=fig.add_axes([fx(26.3),fy(MAP_TOP+MAP_H),MAP_W/W,MAP_H/H])

    base=np.array(Image.open(BASE))

    ax.imshow(base,extent=EXTENT,origin='upper',aspect='auto',interpolation='lanczos',zorder=0)

    a=np.load(_LAYOUT_ROOT / 'results/tracks/absolute_track_point_change.npz')

    change=np.ma.masked_where(np.abs(a['change'])<0.25,a['change'])

    lim=2.0;ALPHA=0.6

    if variant=='overlay':

        m=ax.pcolormesh(a['lon_edges'],a['lat_edges'],change,cmap='RdBu_r',norm=TwoSlopeNorm(0,-lim,lim),shading='flat',alpha=ALPHA,zorder=1,linewidth=0,rasterized=True)

    land=gpd.read_file(PROJECT/'data/raw/GSHHG/GSHHS_shp/i/GSHHS_i_L1.shp',bbox=(EXTENT[0],EXTENT[2],EXTENT[1],EXTENT[3]))

    land.boundary.plot(ax=ax,color='0.40',linewidth=0.3,zorder=2)

    ax.set_xlim(*EXTENT[:2]);ax.set_ylim(*EXTENT[2:])

    ax.set_xticks([100,120,140,160,180]);ax.set_yticks([0,10,20,30,40,50])

    mf.format_geo_axes(ax);ax.tick_params(labelsize=8.6,width=0.8,length=2.4)

    for s in ax.spines.values():s.set_linewidth(0.7)

    lon=lambda x:100+(x-26.3)/4.054;lat=lambda y:(249.1-y)/4.862   # reference PDF points to degrees

    halo=[pe.withStroke(linewidth=2.6,foreground='white',alpha=0.85)]

    ax.add_patch(FancyArrowPatch((lon(280.0),lat(185.2)),(lon(149.3),lat(168.2)),connectionstyle='arc3,rad=0.12',

        arrowstyle='-|>,head_length=4.2,head_width=2.6',color=WEST,lw=3.4,zorder=4,shrinkA=0,shrinkB=0))

    ax.add_patch(FancyArrowPatch((lon(264.0),lat(155.8)),(lon(205.3),lat(66.0)),connectionstyle='arc3,rad=-0.28',

        arrowstyle='-|>,head_length=4.0,head_width=2.4',color=RECURVE,lw=4.4,zorder=4,shrinkA=0,shrinkB=0))

    ax.add_patch(FancyArrowPatch((lon(129.6),lat(128.0)),(lon(129.6),lat(90.0)),arrowstyle='-|>,head_length=4.4,head_width=2.2',

        color=LANDFALL,lw=1.6,zorder=4,shrinkA=0,shrinkB=0))

    t=dict(fontsize=8.6,zorder=5,path_effects=halo)

    ax.text(lon(237.1),lat(212.1),'Low-latitude westward corridor',color=WEST_TEXT,ha='center',va='center',**t)

    ax.text(lon(237.1),lat(227.5),'fewer passages',color=WEST_TEXT,ha='center',va='center',style='italic',**t)

    ax.text(lon(285.7),lat(95.4),'Recurving corridor',color=RECURVE_TEXT,ha='center',va='center',**t)

    ax.text(lon(285.7),lat(110.9),'larger relative share',color=RECURVE_TEXT,ha='center',va='center',style='italic',**t)

    ax.text(lon(136.5),lat(112.2),'poleward shift in',color=LANDFALL,ha='left',va='center',**{**t,'fontsize':8.2})

    ax.text(lon(136.5),lat(122.9),'landfall latitude',color=LANDFALL,ha='left',va='center',**{**t,'fontsize':8.2})

    if variant=='overlay':   # colour bar in the open ocean north of the track-density domain

        ax.add_patch(Rectangle((0.655,0.80),0.338,0.19,transform=ax.transAxes,facecolor='white',edgecolor='none',alpha=0.85,zorder=6))

        cax=ax.inset_axes([0.685,0.845,0.28,0.028],zorder=7)

        cb=fig.colorbar(m,cax=cax,orientation='horizontal',extend='both',ticks=[-2,-1,0,1,2])

        cb.solids.set_alpha(ALPHA);cb.outline.set_linewidth(0.5)

        cax.tick_params(labelsize=7.2,length=1.8,width=0.5,pad=1.2)

        cax.set_title('Track points per year\n1996–2025 minus 1966–1995',fontsize=7.2,pad=2.5,linespacing=1.15)

    ax.text(0.012,0.985,'a',transform=ax.transAxes,ha='left',va='top',fontsize=13,fontweight='bold',zorder=8)



    # ---- (b) flow chart, frame in the reference layout: x 364.7-667.6 pt, y 7.2-247.9 pt (+ legend row)

    bx=fig.add_axes([fx(364.7),fy(276.0),(667.6-364.7)/W,(276.0-7.2)/H]);bx.set_xlim(364.7,667.6);bx.set_ylim(276.0,7.2);bx.set_axis_off()

    bx.add_patch(Rectangle((364.7,7.2),302.9,240.7,facecolor='none',edgecolor=FRAME,lw=0.8,ls=(0,(2.4,1.6))))

    def box(x0,y0,x1,y1,kind,lines,size=8.8):

        style=dict(obs=dict(fc=OBS_FILL,ec=OBS_EDGE,ls='-',tc='white'),assoc=dict(fc=ASSOC_FILL,ec=ASSOC_EDGE,ls='-',tc=ASSOC_TEXT),

                   open=dict(fc='white',ec=OPEN_EDGE,ls=(0,(4.07,1.76)),tc=OPEN_TEXT))[kind]

        bx.add_patch(FancyBboxPatch((x0,y0),x1-x0,y1-y0,boxstyle='round,pad=0,rounding_size=4.5',facecolor=style['fc'],edgecolor=style['ec'],lw=1.1,ls=style['ls'],mutation_aspect=1))

        n=len(lines);step=size*1.22;yc=(y0+y1)/2

        for i,(text,kw) in enumerate(lines):

            bx.text((x0+x1)/2,yc+(i-(n-1)/2)*step,text,ha='center',va='center',color=style['tc'],fontsize=kw.get('size',size),style=kw.get('style','normal'))

    def arrow(p,q,dashed=False,both=False,color=ARROW,lw=1.2):

        bx.add_patch(FancyArrowPatch(p,q,arrowstyle=('<|-|>' if both else '-|>')+',head_length=4.4,head_width=2.2',color=color,lw=lw,

            ls=(0,(3.2,2.0)) if dashed else '-',shrinkA=0,shrinkB=0,zorder=3))

    L=lambda s,**k:(s,k)

    # Row 1: statistical associations

    box(377.9,12.8,654.4,48.0,'assoc',[L('PDO, steering flow, and background-adjusted height'),L('covary with the annual redistribution index (Figure 5)',size=8.0)])

    # Row 2: not resolved by the present diagnostics

    box(368.7,66.0,510.6,113.0,'open',[L('Monthly-mean steering'),L('contribution unresolved'),L('(trajectory-model baseline bias)',size=8.0)])

    box(521.6,66.0,663.6,113.0,'open',[L('Untested processes'),L('submonthly steering regimes,',size=8.0),L('vortex structure, starting position,',size=8.0),L('storm lifetime and termination',size=8.0)])

    # Row 3 and 4: observed changes

    box(390.1,140.0,642.2,174.0,'obs',[L('Relative track composition',size=9.2),L('lower low-latitude share, higher share at 20°–40°N',size=8.0)])

    box(368.7,197.6,510.6,236.0,'obs',[L('Landfall latitude and'),L('coast composition'),L('first landfall +0.9° to +1.4°',size=8.0)])

    box(521.6,197.6,663.6,236.0,'obs',[L('Upper-ocean thermal'),L('exposure along the path'),L('net change unresolved',size=8.0)])

    # Links

    arrow((439.6,48.0),(439.6,65.6),dashed=True);arrow((592.6,48.0),(592.6,65.6),dashed=True)

    arrow((439.6,113.0),(452.0,139.6),dashed=True);arrow((592.6,113.0),(580.0,139.6),dashed=True)

    arrow((516.1,48.4),(516.1,139.6),both=True,color=ASSOC_EDGE,lw=1.3)

    arrow((452.0,174.0),(439.6,197.2));arrow((580.0,174.0),(592.6,197.2))

    # Legend row

    def key(x,kind,label):

        fc,ec,ls={'obs':(OBS_FILL,OBS_EDGE,'-'),'assoc':(ASSOC_FILL,ASSOC_EDGE,'-'),'open':('white',OPEN_EDGE,(0,(2.6,1.2)))}[kind]

        bx.add_patch(FancyBboxPatch((x,256.5),16,9,boxstyle='round,pad=0,rounding_size=2',facecolor=fc,edgecolor=ec,lw=0.9,ls=ls))

        bx.text(x+20,261.0,label,va='center',ha='left',fontsize=8.0,color='#222222')

    key(366.0,'obs','Observed change');key(457.0,'assoc','Statistical association');key(562.0,'open','Unresolved or untested')

    bx.text(366.8,21.5,'b',ha='left',va='center',fontsize=13,fontweight='bold')

    for t in fig.findobj(match=matplotlib.text.Text):t.set_fontfamily(FONT)

    meta={'Title':'Pathway change and levels of evidence','Subject':'Journal of Geophysical Research: Atmospheres Figure 7','Creator':'Matplotlib'}

    stem=out or 'figures/Fig07_outlook_schematic'

    fig.savefig(output_path(stem+'.png'),dpi=600,bbox_inches='tight',pad_inches=0.02)

    fig.savefig(output_path(stem+'.pdf'),bbox_inches='tight',pad_inches=0.02,metadata=meta,dpi=600)

    plt.close(fig)

    if out:print('draft',stem);plt.close(fig);return

    write_json('results/diagnostics/synthesis_figure.json',dict(status='PASS',variant=variant,layout='reference figure (672 x 264 pt) plus a 28-pt legend row',

        shading='results/tracks/absolute_track_point_change.npz, cells with |change|<0.15 masked, colour limits +/-2 track points per year (extend both)',

        evidence_levels=dict(observed=['Relative track composition','Landfall latitude and coast composition','Upper-ocean thermal exposure along the path'],

            association=['PDO, steering flow, and background-adjusted height covary with the annual redistribution index'],

            unresolved_or_untested=['Monthly-mean steering contribution','Submonthly steering regimes, vortex structure, starting position, storm lifetime and termination']),

        removed=['Reconfiguration of subtropical-circulation geometry','PDO-to-circulation causal arrow','Displacement of the mean steering flow (0.14 sigma zonal, 0.08 sigma meridional)']))

    print('Figure 7 written')



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')



    {'download':download,'build':build,'draw':draw}[sys.argv[1]](*sys.argv[2:])
