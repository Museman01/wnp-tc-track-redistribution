'Synthesis figure (Figure 7): panel b adds the circulation association (double-headed link) and untested processes (dotted box) to the partitions, interpretation, and expressions.'
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import sys,time

import numpy as np

from PIL import Image

import analysis.figures.relief_basemap as mb

from core.paths import PROJECT,output_path



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



def build(output):

    out=(PROJECT / output).resolve()
    if not out.is_relative_to(PROJECT.resolve()): raise ValueError('Output must remain within the project')
    if out.exists(): raise FileExistsError('Refusing to overwrite an existing basemap')
    out.parent.mkdir(parents=True,exist_ok=True)
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

    Image.fromarray(grey).save(out)

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
    from pathlib import Path
    if out is None or not (PROJECT / out).resolve().is_relative_to(PROJECT.resolve()):
        raise ValueError('Provide an output stem within the project.')
    if any(Path(str((PROJECT / out).resolve()) + suffix).exists() for suffix in ('.png', '.pdf')):
        raise FileExistsError('Refusing to overwrite an existing figure.')

    import matplotlib

    matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    import matplotlib.patheffects as pe

    from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Rectangle

    from matplotlib.colors import TwoSlopeNorm

    import geopandas as gpd

    import analysis.figures.main_figures as mf   # font set-up and geographic tick format

    FONT=mf.FONT

    W,H=672.0,428.0                     # reference figure width 672 pt; height for the colour bar below panel a and two legend rows in panel b

    fig=plt.figure(figsize=(W/72,H/72))

    fx=lambda x:x/W;fy=lambda y:1-y/H  # PDF points (origin top-left) to figure fraction



    # ---- (a) map, same frame in the reference layout: x 26.3-350.6 pt, y 15.75-239.35 pt; colour bar below the map
    DY=(H-320.0)/2                    # centre panel a against the taller panel b
    ax=fig.add_axes([fx(26.3),fy(MAP_TOP+DY+MAP_H),MAP_W/W,MAP_H/H])
    base=np.array(Image.open(BASE))
    ax.imshow(base,extent=EXTENT,origin='upper',aspect='auto',interpolation='lanczos',zorder=0)
    a=np.load(_LAYOUT_ROOT / 'results/tracks/absolute_track_point_change.npz')
    change=np.ma.masked_where(np.abs(a['change'])<0.25,a['change'])
    lim=2.0;ALPHA=0.6
    m=ax.pcolormesh(a['lon_edges'],a['lat_edges'],change,cmap='RdBu_r',norm=TwoSlopeNorm(0,-lim,lim),shading='flat',alpha=ALPHA,zorder=1,linewidth=0,rasterized=True)
    land=gpd.read_file(PROJECT/'data/raw/GSHHG/GSHHS_shp/i/GSHHS_i_L1.shp',bbox=(EXTENT[0],EXTENT[2],EXTENT[1],EXTENT[3]))
    land.boundary.plot(ax=ax,color='0.40',linewidth=0.3,zorder=2)
    ax.set_xlim(*EXTENT[:2]);ax.set_ylim(*EXTENT[2:])
    ax.set_xticks([100,120,140,160,180]);ax.set_yticks([0,10,20,30,40,50])
    mf.format_geo_axes(ax);ax.tick_params(labelsize=9.0,width=0.8,length=2.4)
    for s in ax.spines.values():s.set_linewidth(0.7)
    lon=lambda x:100+(x-26.3)/4.054;lat=lambda y:(249.1-y)/4.862   # reference PDF points to degrees
    halo=[pe.withStroke(linewidth=2.6,foreground='white',alpha=0.85)]
    ax.add_patch(FancyArrowPatch((lon(280.0),lat(185.2)),(lon(149.3),lat(168.2)),connectionstyle='arc3,rad=0.12',
        arrowstyle='-|>,head_length=4.2,head_width=2.6',color=WEST,lw=3.4,zorder=4,shrinkA=0,shrinkB=0))
    ax.add_patch(FancyArrowPatch((lon(264.0),lat(155.8)),(lon(205.3),lat(66.0)),connectionstyle='arc3,rad=-0.28',
        arrowstyle='-|>,head_length=4.0,head_width=2.4',color=RECURVE,lw=4.4,zorder=4,shrinkA=0,shrinkB=0))
    ax.add_patch(FancyArrowPatch((lon(129.6),lat(128.0)),(lon(129.6),lat(90.0)),arrowstyle='-|>,head_length=4.4,head_width=2.2',
        color=LANDFALL,lw=1.6,zorder=4,shrinkA=0,shrinkB=0))
    t=dict(fontsize=9.0,zorder=5,path_effects=halo)
    ax.text(lon(237.1),lat(212.1),'Low-latitude westward corridor',color=WEST_TEXT,ha='center',va='center',**t)
    ax.text(lon(237.1),lat(227.5),'fewer passages',color=WEST_TEXT,ha='center',va='center',style='italic',**t)
    ax.text(lon(285.7),lat(95.4),'Recurving corridor',color=RECURVE_TEXT,ha='center',va='center',**t)
    ax.text(lon(285.7),lat(110.9),'larger relative share',color=RECURVE_TEXT,ha='center',va='center',style='italic',**t)
    ax.text(lon(136.5),lat(112.2),'Poleward shift in',color=LANDFALL,ha='left',va='center',**{**t,'fontsize':8.6})
    ax.text(lon(136.5),lat(122.9),'landfall latitude',color=LANDFALL,ha='left',va='center',**{**t,'fontsize':8.6})
    cax=fig.add_axes([fx(98.0),fy(266.0+DY),180.0/W,6.5/H])
    cb=fig.colorbar(m,cax=cax,orientation='horizontal',extend='both',ticks=[-2,-1,0,1,2])
    cb.solids.set_alpha(ALPHA);cb.outline.set_linewidth(0.5)
    cax.tick_params(labelsize=9.0,length=1.8,width=0.5,pad=1.2)
    cb.ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,_:f'{v:g}'.replace('-','−')))
    cax.set_xlabel('Difference in mean annual track-point count\nper 2.5° × 2.5° grid cell\n1996–2025 minus 1966–1995',fontsize=9.0,labelpad=2.5,linespacing=1.2)
    ax.text(0.012,0.985,'a',transform=ax.transAxes,ha='left',va='top',fontsize=13,fontweight='bold',zorder=8)

    # ---- (b) shared result, statistical analyses, interpretation, untested processes, and expressions
    BS,TS=9.8,10.4
    SB,ST=9.4,10.0                      # body and title sizes in the three narrow boxes
    YB=H-3.0
    bx=fig.add_axes([fx(364.7),fy(YB),(667.6-364.7)/W,(YB-4.0)/H])
    bx.set_xlim(364.7,667.6);bx.set_ylim(YB,4.0);bx.set_axis_off()
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
        for text,offset,size,weight in items:
            bx.text((x0+x1)/2,y0+offset,text,ha='center',va='center',
                fontsize=size,fontweight=weight,color=tc)
    def segment(x0,y0,x1,y1,ls='solid'):
        bx.plot([x0,x1],[y0,y1],color='#62696c',lw=.85,linestyle=ls,solid_capstyle='butt',zorder=1)
    def arrow(x0,y0,x1,y1,style='-|>',ls='solid'):
        bx.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle=style,
            mutation_scale=8,lw=.85,color='#62696c',linestyle=ls,shrinkA=0,shrinkB=0,zorder=2))
    X0,X1=368.7,663.6
    GAP=6.0;CW=(X1-X0-2*GAP)/3
    c1=(X0,X0+CW);c2=(X0+CW+GAP,X0+2*CW+GAP);c3=(X1-CW,X1)
    m1=sum(c1)/2;m2=sum(c2)/2;m3=sum(c3)/2;m12=(c1[0]+c2[1])/2
    L,R=440.45,591.85                  # centres of the two expression boxes, as in the previous layout
    # row 1: shared observed result
    box(X0,12,X1,57,'observed',[
        ('Shared shift toward higher-latitude pathways',16,10.8,'bold'),
        ('Relative track composition across four agencies',32,BS,'normal')])
    # row 2: two partitions (arrows) and one association (double-headed line, no causal direction)
    Y2,Y2b=92,214
    segment(m12,57,m12,66);segment(m1,66,m2,66)
    arrow(m1,66,m1,Y2-2);arrow(m2,66,m2,Y2-2)
    bx.text(m12,78,'Partitions',ha='center',va='center',fontsize=BS,color='#555d61')
    arrow(m3,57,m3,Y2-2,style='<|-|>')
    bx.text(m3-6,74,'Association',ha='right',va='center',fontsize=BS,color='#555d61')
    box(c1[0],Y2,c1[1],Y2b,'partition',[
        ('Genesis-bin',13,ST,'bold'),('partition',26,ST,'bold'),
        ('74%–80%',50,16,'bold'),
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
    Y3,Y3b=252,320
    segment(m1,Y2b,m1,Y2b+8);segment(m2,Y2b,m2,Y2b+8);segment(m1,Y2b+8,m2,Y2b+8)
    arrow(m12,Y2b+8,m12,Y3-2)
    bx.text(m12+6,Y2b+22,'Interpretation',ha='left',va='center',fontsize=BS,color='#555d61')
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
    Y4,Y4b=356,420
    segment(m12,Y3b,m12,Y3b+10);segment(L,Y3b+10,R,Y3b+10)
    arrow(L,Y3b+10,L,Y4-2);arrow(R,Y3b+10,R,Y4-2)
    bx.text((L+R)/2,Y3b+23,'Along-track expressions',ha='center',va='center',fontsize=BS,color='#555d61')
    box(368.7,Y4,512.2,Y4b,'expression',[
        ('Coastal displacement',15,TS,'bold'),
        ('Higher mean first-landfall',34,BS,'normal'),
        ('latitude in the later period',49,BS,'normal')])
    box(520.1,Y4,663.6,Y4b,'expression',[
        ('Upper-ocean exposure',12,TS,'bold'),
        ('Path and ocean changes',27,BS,'normal'),
        ('act in opposite directions;',41,BS,'normal'),
        ('net change unresolved',55,BS,'normal')])
    bx.text(355.8,6,'b',ha='left',va='top',fontsize=13,fontweight='bold')
    for t in fig.findobj(match=matplotlib.text.Text):t.set_fontfamily(FONT)

    meta={'Title':'Pathway change, statistical analyses, and untested processes','Subject':'Journal of Geophysical Research: Atmospheres Figure 7','Creator':'Matplotlib'}

    stem=str(out)

    fig.savefig(output_path(stem+'.png'),dpi=600,bbox_inches='tight',pad_inches=0.02)
    fig.savefig(output_path(stem+'.pdf'),bbox_inches='tight',pad_inches=0.02,metadata=meta,dpi=600)

    plt.close(fig)
    print('Figure saved:',stem)

if __name__=='__main__':
    if '--execute' not in sys.argv:
        raise SystemExit('Explicit execution required: add --execute.')
    sys.argv.remove('--execute')
    {'download':download,'build':build,'draw':draw}[sys.argv[1]](*sys.argv[2:])
