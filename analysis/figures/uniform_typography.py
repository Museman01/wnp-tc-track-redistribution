"""Uniform figure lettering on the original scientific panel layouts."""
from pathlib import Path
import argparse
import copy
from collections import Counter
import hashlib
import importlib.util
import json
import re
import sys

from core.paths import PROJECT_ROOT as PROJECT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text

OUTPUT = None
CONFIG = json.loads((PROJECT/'data/figure_inputs/typography/display_parameters.json').read_text())

def route(relative):
    path = OUTPUT / Path(relative).name
    return path

def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, PROJECT/'analysis/figures'/filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    # Route all scientific renderers to the explicitly selected figure directory.
    mod.output_path = lambda relative: route(relative) if str(relative).startswith('figures/') else PROJECT/relative
    return mod


def array_hash(value):
    arr = np.ma.asarray(value)
    if arr.dtype.kind in 'OUS':
        raw = repr(arr.tolist()).encode()
    else:
        raw = np.ascontiguousarray(arr.filled(np.nan) if np.ma.isMaskedArray(arr) and np.any(arr.mask) else arr).tobytes()
    return [list(arr.shape), str(arr.dtype), hashlib.sha256(raw).hexdigest()]


def scientific_payload(fig):
    result = []
    for ax in fig.axes:
        record = dict(xlim=list(ax.get_xlim()), ylim=list(ax.get_ylim()), lines=[], patches=[], images=[], collections=[])
        for line in ax.lines:
            record['lines'].append([array_hash(line.get_xdata()), array_hash(line.get_ydata()),
                                    line.get_color(), line.get_linestyle(), line.get_marker()])
        for patch in ax.patches:
            record['patches'].append([array_hash(patch.get_path().vertices),
                                      list(patch.get_facecolor()), list(patch.get_edgecolor())])
        for image in ax.images:
            record['images'].append([array_hash(image.get_array()), list(image.get_extent()),
                                     image.get_cmap().name, list(image.get_clim())])
        for collection in ax.collections:
            record['collections'].append([array_hash(collection.get_offsets()),
                 array_hash(collection.get_array()) if collection.get_array() is not None else None,
                 [array_hash(p.vertices) for p in collection.get_paths()],
                 array_hash(collection.get_facecolors()), array_hash(collection.get_edgecolors()),
                 list(collection.get_clim())])
        # Axis tick selection/positions are presentation; their ranges are fixed above.
        words = [ax.get_xlabel(), ax.get_ylabel(), ax.get_title()] + [t.get_text() for t in ax.texts]
        legend = ax.get_legend()
        if legend:
            words += [t.get_text() for t in legend.get_texts()]
        record['text'] = sorted(Counter(re.sub(r'\s+', '', s) for s in words if s).items())
        result.append(record)
    return result



PROFILES={}

def save_uniform(fig,stem,scale=1,**kwargs):
    if scale != 1:sys.modules['analysis.figures.main_figures'].scale_typography(fig,scale)
    fig.canvas.draw()
    before=scientific_payload(fig)
    positions=[list(ax.get_position().bounds) for ax in fig.axes]
    sizes=list(fig.get_size_inches())
    fig.set_layout_engine(None)
    from matplotlib.ticker import FixedLocator
    for ax,pos in zip(fig.axes,positions):
        ax.set_position(pos)
        ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        ax.yaxis.set_major_locator(FixedLocator(ax.get_yticks()))
    renderer=fig.canvas.get_renderer()
    legends=list(fig.legends)+[ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]
    for ax in fig.axes:
        leg=ax.get_legend()
        if leg and leg._loc==0:
            b=leg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
            leg.set_loc('upper right');leg.set_bbox_to_anchor((b.x1,b.y1),transform=ax.transAxes);leg.borderaxespad=0
    base=13.
    if stem.startswith('Fig02'):base=14.
    if stem.startswith(('Fig03','Fig04')):base=12.5
    if stem.startswith('Fig05'):base=14.
    if stem.startswith('Fig07'):base=8.6
    if stem.startswith('Fig02'):
        for ax in fig.axes:
            for t in ax.texts:
                if t.get_text()=='p=0.0015':t.set_transform(t.get_transform()+matplotlib.transforms.ScaledTranslation(-4/72,0,fig.dpi_scale_trans))
    if stem.startswith('Fig03'):
        ax=next(ax for ax in fig.axes if ax.get_legend() and 'Primary' in [t.get_text() for t in ax.get_legend().get_texts()]);leg=ax.get_legend();anchor=leg.get_bbox_to_anchor().transformed(ax.transAxes.inverted())
        leg.set_bbox_to_anchor((anchor.x0-.04,anchor.y0-.025),transform=ax.transAxes)
    if stem.startswith('FigS11'):
        ax=fig.axes[3];leg=ax.get_legend();anchor=leg.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transAxes.inverted())
        leg.set_loc('upper right');leg.set_bbox_to_anchor((anchor.x1,anchor.y1-.12),transform=ax.transAxes);leg.borderaxespad=0
    if stem.startswith('FigS09'):base=12.5
    if stem.startswith('FigS15'):base=11.5
    if stem.startswith('FigS19'):base=11.5
    targets={role:base for role in ['body','legend','annotation','panel','symbol','colorbar']}
    removed=[]
    for i,ax in enumerate(fig.axes):
        for loc in ([] if hasattr(ax,'_colorbar') else ['left','center','right']):
            value=ax.get_title(loc=loc)
            if value:removed.append(dict(axis=i,location=loc,text=value));ax.set_title('',loc=loc)
    if fig._suptitle is not None:removed.append(dict(axis=None,text=fig._suptitle.get_text()));fig._suptitle.set_text('')
    for i,ax in enumerate(fig.axes):
        for text in list(ax.texts):
            if text.get_text() in ['Total change','Genesis distribution','Conditional track field','Raw','Detrended']:
                removed.append(dict(axis=i,location='panel heading',text=text.get_text()));text.remove()
    expected=scientific_payload(fig)
    from matplotlib.quiver import QuiverKey
    keys=[v for ax in fig.axes for v in ax.get_children() if isinstance(v,QuiverKey)]
    if stem.startswith(('FigS16','FigS17')):
        for key in keys:key.X-=.025
    all_text=list(fig.findobj(match=Text))+[v.text for v in keys]
    role_items=[]
    for text in all_text:
        if text.get_visible() and text.get_text().strip():
            value=text.get_text();role='body'
            if any(text in leg.get_texts() for leg in legends):role='legend'
            elif re.fullmatch(r'\(?[a-z]\)?',value) and any(text in ax.texts for ax in fig.axes):role='panel'
            elif value=='*':role='symbol'
            elif any(text in ax.texts for ax in fig.axes):role='annotation'
            if stem=='Fig07_outlook_schematic' and any(text in child.findobj(match=Text) for ax in fig.axes for child in ax.child_axes):role='colorbar'
            text.set_fontsize(targets[role]);role_items.append([value,role,targets[role]])
    if stem.startswith('Fig02'):
        for ax in fig.axes:
            for t in ax.texts:
                if t.get_text()=='p=0.0015':t.set_transform(t.get_transform()+matplotlib.transforms.ScaledTranslation(-4/72,0,fig.dpi_scale_trans))
    if stem.startswith('Fig03'):
        ax=next(ax for ax in fig.axes if ax.get_legend() and 'Primary' in [t.get_text() for t in ax.get_legend().get_texts()]);leg=ax.get_legend();anchor=leg.get_bbox_to_anchor().transformed(ax.transAxes.inverted())
        leg.set_bbox_to_anchor((anchor.x0-.04,anchor.y0-.025),transform=ax.transAxes)
    if stem.startswith('FigS11'):
        ax=fig.axes[3];leg=ax.get_legend();anchor=leg.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transAxes.inverted())
        leg.set_loc('upper right');leg.set_bbox_to_anchor((anchor.x1,anchor.y1-.12),transform=ax.transAxes);leg.borderaxespad=0
    if stem.startswith('FigS09'):
        for ax in fig.axes:
            leg=ax.get_legend()
            if leg:
                for t in leg.get_texts():
                    if t.get_text() in ['Track-density term','Regional-speed term','Covariance term']:
                        t.set_text(t.get_text().replace(' term','\nterm'))
                leg.columnspacing=.6
    if stem.startswith('Fig07'):
        for ax in fig.axes:
            for child in ax.child_axes:
                if child.get_title():child.title.set_x(.34)
    # Uniform panel letters can require a few points of clearance at the same corner.
    if stem.startswith('FigS14'):
        for ax in fig.axes:
            for t in ax.texts:
                if re.fullmatch('[a-c]',t.get_text()):t.set_x(.005)
    if stem.startswith('FigS12'):
        for ax in fig.axes:
            for t in ax.texts:
                if re.fullmatch('[a-h]',t.get_text()):t.set_position((.012,1.045));t.set_clip_on(False)
    fig.canvas.draw()
    assert expected==scientific_payload(fig),(stem,'unintended scientific change')
    assert np.allclose(positions,[ax.get_position().bounds for ax in fig.axes])
    # A small tick-label gap avoids southwest-corner collisions without moving panels.
    for ax in fig.axes:
        labels=[v.get_text() for v in ax.get_xticklabels()]
        if any('°' in v for v in labels):
            pad=ax.xaxis.get_major_ticks()[0].get_pad() if ax.xaxis.get_major_ticks() else 3.5
            ax.tick_params(axis='x',pad=pad+2)
    if stem.startswith('FigS15'):
        for ax in fig.axes[:3]:
            ax.tick_params(axis='x',pad=3.5)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();inactive=set()
    for ax in fig.axes:
        for axis in [ax.xaxis,ax.yaxis]:
            if not ax.axison:
                inactive.update([axis.label,axis.get_offset_text()])
                for tick in axis.get_major_ticks()+axis.get_minor_ticks():inactive.update([tick.label1,tick.label2])
            active=set(axis._update_ticks())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if tick not in active:inactive.update([tick.label1,tick.label2])
    boxes=[]
    for text in all_text:
        if text not in inactive and text.get_visible() and text.get_text().strip():boxes.append((text.get_text(),text.get_window_extent(renderer)))
    conflicts=[]
    for i,(a,ba) in enumerate(boxes):
        for b,bb in boxes[i+1:]:
            hit=matplotlib.transforms.Bbox.intersection(ba,bb)
            if hit is not None and hit.width>.5 and hit.height>.5:conflicts.append([a,b,[round(v,2) for v in hit.bounds]])
    fixed=stem.startswith(('FigS15','FigS16','FigS17','FigS18'))
    crop=fig.bbox_inches if fixed else matplotlib.transforms.Bbox.union([fig.get_tightbbox(renderer)]+[box.transformed(fig.dpi_scale_trans.inverted()) for _,box in boxes]).padded(.16 if stem.startswith('Fig04') else .10)
    options={} if fixed else dict(bbox_inches=crop)
    export_dpi=CONFIG[stem]['dpi']
    ORIGINAL_SAVEFIG(fig,route('figures/'+stem+'.png'),dpi=export_dpi,**options)
    ORIGINAL_SAVEFIG(fig,route('figures/'+stem+'.pdf'),dpi=export_dpi,**options,metadata={'Title':stem.replace('_',' '),'Author':'Daizhao Huang; Yuanfu Yue'})
    plt.close(fig)
    print(stem, flush=True)

from matplotlib.figure import Figure
ORIGINAL_SAVEFIG=Figure.savefig

class FinishedFigure(Exception):pass

def main():
    global OUTPUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('figure',choices=['2','3','4','5','6','7','S01','S02','S03','S06','S07','S09','S11','S12','S13','S14','S15','S16','S17','S18','S19'])
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--execute', action='store_true')
    a=parser.parse_args()
    if not a.execute:parser.error('Use --execute to draw figures')
    OUTPUT=a.output_dir.resolve()
    if not OUTPUT.is_relative_to(PROJECT.resolve()):parser.error('Output directory must remain within the project')
    OUTPUT.mkdir(parents=True,exist_ok=True)
    # Each selected renderer declares its output names; reject an occupied directory.
    if any(OUTPUT.iterdir()):parser.error('Select an empty output directory')
    mf=module('analysis.figures.main_figures','main_figures.py')
    if a.figure in ['2','3','4','5','6']:
        mf.save=save_uniform;getattr(mf,'figure'+a.figure)()
    elif a.figure=='7':
        mod=module('analysis.figures.synthesis_figure','synthesis_figure.py')
        def intercept(fig,*args,**kw):save_uniform(fig,'Fig07_outlook_schematic');raise FinishedFigure()
        Figure.savefig=intercept
        try:mod.draw(out=str(OUTPUT/'Fig07_outlook_schematic'))
        except FinishedFigure:pass
        finally:Figure.savefig=ORIGINAL_SAVEFIG
    elif a.figure in ['S14','S15','S16','S17','S18','S19']:
        d=module('analysis.figures.diagnostic_figures','diagnostic_figures.py')
        targets={'S14':('absolute','FigS14_absolute_passage_recurvature'),
                 'S15':('bam','FigS15_bam_baseline_skill'),
                 'S16':('steering','FigS16_steering_layer_regressions'),
                 'S17':('steering','FigS17_circulation_geometry_diagnostics'),
                 'S18':('hko','FigS18_hko_common_window'),
                 'S19':('d26','FigS19_d26_parallel_exposure')}
        renderer,stem=targets[a.figure]
        d.save=lambda fig,name:save_uniform(fig,name) if name==stem else plt.close(fig)
        getattr(d,renderer)()
    else:
        si=module('analysis.figures.si_four_agency_figures','si_four_agency_figures.py')
        si.save=lambda fig,stem,key,title,scale:save_uniform(fig,stem,scale=scale)
        getattr(si,a.figure)()

if __name__=='__main__':main()
