"""Supporting Information Figure S5 redrawn with the map at true scale at 20°N.



The plotting code of the submitted Figure S5 is not part of the project; its values are

therefore read from the submitted vector PDF itself and re-plotted unchanged:

  * panel a: the 32 x 16 grid of 2.5-degree cells is stored as filled rectangles whose

    viridis colours (vmin 0, vmax 25) map back to integer storm counts;

  * panel b: marker centres give mean LMI (kt) per 5-degree genesis-latitude bin, and the

    n labels give the storm numbers.

Extracted values: results/tracks/figS05_lmi_location_values.csv and .._genesis_bins.csv."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import re

import numpy as np,pandas as pd

import pymupdf

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib import colormaps

from core.paths import output_path,write_json,PROJECT

import analysis.figures.main_figures as mf



SRC=_LAYOUT_ROOT / 'data/figure_inputs/FigS05_lmi_location_reference.pdf'

LAND,EDGE='#ECE9E1','#606569'



def extract():

    page=pymupdf.open(SRC)[0];lut=colormaps['viridis'](np.arange(256))[:,:3]

    cells=[]

    for dr in page.get_drawings():

        r=dr['rect'];f=dr.get('fill')

        if f and dr['type']=='f' and 11<r.width<13.5 and 14<r.height<16.5 and r.x0>=35 and r.x1<=430:

            k=int(np.argmin(((lut-np.array(f))**2).sum(1)));cells.append((r.x0,r.y1,k))

    assert len(cells)==512

    x0=min(c[0] for c in cells);y1=max(c[1] for c in cells);w=(max(c[0] for c in cells)-x0)/31;h=(y1-min(c[1] for c in cells))/15

    grid=np.full((16,32),np.nan)

    for x,y,k in cells:

        i=int(round((y1-y)/h));j=int(round((x-x0)/w));grid[i,j]=k

    # Matplotlib maps a value v to LUT index floor(256*v/25), capped at 255.

    lookup={min(int(np.floor(256*n/25)),255):n for n in range(26)}

    assert not np.isnan(grid).any() and set(np.unique(grid).astype(int))<=set(lookup)

    counts=np.vectorize(lambda k:lookup[int(k)])(grid).astype(float)

    # panel b: axis transforms from tick marks (50/75/100 kt; 0/40 deg)

    ticks=sorted({round(dr['rect'].x0,3) for dr in page.get_drawings() if dr['type']=='fs' and dr['rect'].width<0.01 and dr['rect'].x0>440 and abs(dr['rect'].y1-34.88)<0.1})

    assert len(ticks)==3;kt_per_pt=25/(ticks[1]-ticks[0])

    marks=[dr['rect'] for dr in page.get_drawings() if dr['type']=='fs' and dr['rect'].x0>440 and 4<dr['rect'].width<5]

    labels=[(s['bbox'][1],int(s['text'][2:])) for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if re.fullmatch(r'n=\d+',s['text'])]

    yaxis=(278.24,34.88)

    rows=[]

    for r in marks:

        cx=(r.x0+r.x1)/2;cy=(r.y0+r.y1)/2

        lat=(yaxis[0]-cy)/(yaxis[0]-yaxis[1])*40;kt=50+(cx-ticks[0])*kt_per_pt

        n=min(labels,key=lambda t:abs(t[0]+6-cy))[1]

        rows.append(dict(bin_center_deg=round(lat,2),mean_lmi_kt=round(kt,2),n_storms=n))

    bins=pd.DataFrame(rows).sort_values('bin_center_deg')

    xlim=(50-(ticks[0]-469.65)*kt_per_pt,100+(597.8-ticks[2])*kt_per_pt)

    return counts,bins,xlim



def main():

    counts,bins,xlim=extract()

    lon=np.arange(100,182.5,2.5);lat=np.arange(0,42.5,2.5)

    pd.DataFrame(counts,index=[f'{a:g}-{a+2.5:g}N' for a in lat[:-1]],columns=[f'{a:g}-{a+2.5:g}E' for a in lon[:-1]]).to_csv(_LAYOUT_ROOT / 'results/tracks/figS05_lmi_location_values.csv')

    bins.to_csv(_LAYOUT_ROOT / 'results/tracks/figS05_lmi_genesis_bins.csv',index=False)

    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':[mf.FONT,'DejaVu Sans'],'pdf.fonttype':42})

    W=627.7;H=300.0

    fig=plt.figure(figsize=(W/72,H/72))

    mx0,mx1=35.17,429.49;mh=(mx1-mx0)*0.5/np.cos(np.deg2rad(20));mtop=8.0

    ax=fig.add_axes([mx0/W,1-(mtop+mh)/H,(mx1-mx0)/W,mh/H])

    m=ax.pcolormesh(lon,lat,counts,cmap='viridis',vmin=0,vmax=25,shading='flat',zorder=1,rasterized=False)

    import geopandas as gpd

    land=gpd.read_file(PROJECT/'data/raw/GSHHG/GSHHS_shp/i/GSHHS_i_L1.shp',bbox=(100,0,180,40))

    land.plot(ax=ax,facecolor=LAND,edgecolor=EDGE,linewidth=0.35,zorder=3)

    ax.set_xlim(100,180);ax.set_ylim(0,40);ax.set_aspect(1/np.cos(np.deg2rad(20)),adjustable='box')

    ax.set_xticks([100,120,140,160,180]);ax.set_yticks([0,10,20,30,40]);mf.format_geo_axes(ax);ax.tick_params(labelsize=9.4,width=0.8)

    for s in ['top','right']:ax.spines[s].set_visible(True)

    ax.text(0.012,0.975,'a',transform=ax.transAxes,ha='left',va='top',fontsize=11.5,fontweight='bold',zorder=6)

    cax=fig.add_axes([92.88/W,1-(mtop+mh+44)/H,(372.24-92.88)/W,10.8/H])

    cb=fig.colorbar(m,cax=cax,orientation='horizontal');cb.set_label('LMI-location density',fontsize=10.8);cax.tick_params(labelsize=9.4)

    bx=fig.add_axes([469.65/W,1-(mtop+mh)/H,(597.8-469.65)/W,mh/H])

    bx.plot(bins.mean_lmi_kt,bins.bin_center_deg,color='#4C72B0',lw=1.45,marker='o',ms=4.5)

    for _,r in bins.iterrows():bx.text(r.mean_lmi_kt+2.2,r.bin_center_deg,f'n={int(r.n_storms)}',va='center',fontsize=10)

    bx.set_xlim(*xlim);bx.set_ylim(0,40);bx.set_yticks(range(0,45,5));bx.set_xticks([50,75,100])

    bx.xaxis.tick_top();bx.xaxis.set_label_position('top');bx.set_xlabel('Mean LMI (kt)',fontsize=10.8);bx.set_ylabel('Genesis-latitude bin (°N)',fontsize=10.8)

    for s in ['right','bottom']:bx.spines[s].set_visible(False)

    bx.tick_params(labelsize=9.4,width=0.8)

    bx.text(0.97,0.975,'b',transform=bx.transAxes,ha='right',va='top',fontsize=11.5,fontweight='bold')

    for t in fig.findobj(match=matplotlib.text.Text):t.set_fontfamily(mf.FONT)

    stem='figures/FigS05_lmi_location'

    fig.savefig(output_path(stem+'.png'),dpi=600,bbox_inches='tight',pad_inches=0.08)

    fig.savefig(output_path(stem+'.pdf'),bbox_inches='tight',pad_inches=0.08,metadata={'Creator':'Matplotlib'})

    plt.close(fig)

    write_json('results/diagnostics/figS05_redraw.json',dict(status='PASS',source='values read from the submitted Figure S5 vector PDF',cells=int(counts.size),

        total_count=int(counts.sum()),max_count=int(counts.max()),bins=bins.to_dict('records'),map_aspect='1/cos(20 deg)'))

    print('S5',int(counts.sum()),bins.to_dict('records'))



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
