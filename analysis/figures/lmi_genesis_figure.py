"""Plot LMI location and genesis-latitude statistics from structured inputs."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT


import numpy as np,pandas as pd


import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt


from core.paths import output_path,PROJECT

import analysis.figures.main_figures as mf




LAND,EDGE='#ECE9E1','#606569'



def extract():
    counts = pd.read_csv(PROJECT / 'results/tracks/figS05_lmi_location_values.csv', index_col=0).to_numpy(dtype=float)
    bins = pd.read_csv(PROJECT / 'results/tracks/figS05_lmi_genesis_bins.csv')
    return counts, bins, (30.959451583358803, 119.87594016263766)


def main(stem):
    if not (PROJECT / stem).resolve().is_relative_to(PROJECT.resolve()):
        raise ValueError('Provide an output stem within the project')
    if any((PROJECT / (stem + suffix)).exists() for suffix in ('.png', '.pdf')):
        raise FileExistsError('Refusing to overwrite an existing figure')
    counts,bins,xlim=extract()
    lon=np.arange(100,182.5,2.5);lat=np.arange(0,42.5,2.5)

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

    fig.savefig(output_path(stem+'.png'),dpi=600,bbox_inches='tight',pad_inches=0.08)

    fig.savefig(output_path(stem+'.pdf'),bbox_inches='tight',pad_inches=0.08,metadata={'Creator':'Matplotlib'})

    plt.close(fig)

    print('S5',int(counts.sum()),bins.to_dict('records'))



if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--output',required=True,help='Output stem within the project')
    args=parser.parse_args()
    if not args.execute:parser.error('add --execute to draw the figure')
    main(args.output)
