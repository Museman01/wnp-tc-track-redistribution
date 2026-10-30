'Study-region map using a reference vector PDF, GSHHG coastline and Maps-For-Free relief. Terrain tint and saturation are both 0.30; PNG output is rendered at 600 dpi.'

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import io

from pathlib import Path

import numpy as np

import pymupdf

from PIL import Image

from core.paths import output_path,write_json



SRC=_LAYOUT_ROOT / 'data/figure_inputs/study_region_reference.pdf'

FONT=next((str(p) for p in [Path('C:/Windows/Fonts/arial.ttf'), Path('/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf')] if p.is_file()), None)

if FONT is None: raise FileNotFoundError('Arial or Liberation Sans is required for the existing latitude labels')

LINE_COLOR=(0.2,0.2,0.2)



def recolour(rgb,sigma=6,post=1.2):

    """Keep only terrain shading: luminance divided by its local mean isolates the

    hill-shade texture from the Maps-For-Free elevation/depth tint."""

    from PIL import ImageFilter

    a=rgb.astype(float);r,g,b=a[...,0],a[...,1],a[...,2]

    sea=(b-np.maximum(r,g))>18

    L=0.299*r+0.587*g+0.114*b

    smooth=np.array(Image.fromarray(L.clip(0,255).astype(np.uint8)).filter(ImageFilter.GaussianBlur(sigma)),float)+1

    shade=L/smooth

    # Light smoothing of the shade suppresses single-pixel texture from the source tiles.

    if post:shade=np.array(Image.fromarray(((shade/2).clip(0,1)*255).astype(np.uint8)).filter(ImageFilter.GaussianBlur(post)),float)/255*2

    out=np.zeros_like(a)

    for mask,base,gain in [(~sea,np.array([212,211,207]),0.95),(sea,np.array([226,234,240]),0.5)]:

        k=1+gain*(shade[mask]-1)

        out[mask]=base*np.clip(k,0.55,1.18)[:,None]

    return out.clip(0,255).astype(np.uint8),sea.mean()



TERRAIN_TINT=0.30

TERRAIN_SATURATION=0.30



def terrain_tint(rgb,strength=TERRAIN_TINT,saturation=TERRAIN_SATURATION):

    """Reduce original basemap colour saturation, then blend with white.

    Both factors are relative to the unmodified Maps-For-Free source raster."""

    from PIL import ImageEnhance

    if not (0 <= strength <= 1 and 0 <= saturation <= 1):

        raise ValueError('Terrain opacity and saturation must be between zero and one')

    original=rgb.astype(float)

    sea=(original[...,2]-np.maximum(original[...,0],original[...,1]))>18

    toned=np.asarray(ImageEnhance.Color(Image.fromarray(rgb)).enhance(saturation),dtype=float)

    out=np.rint(255.0+strength*(toned-255.0)).clip(0,255).astype(np.uint8)

    return out,sea.mean()



RIVER_RG=b'0.0431372549 0.4705882353 0.7058823529 RG'



def strip_rivers(doc,page):

    """Remove the river line groups (one q...Q block per river width class) from the page stream."""

    xrefs=page.get_contents();assert len(xrefs)==1

    data=doc.xref_stream(xrefs[0]);removed=0

    while True:

        i=data.find(RIVER_RG)

        if i<0:break

        start=data.rfind(b'Q q',0,i)+2;end=data.find(b'\nQ q',i)

        assert start>1 and end>i and data[start:start+1]==b'q'

        data=data[:start]+data[end+3:];removed+=1   # keeps one "Q q" joint between neighbours

    doc.update_stream(xrefs[0],data);return removed



def gaps(page,y,pad=2.5):

    """Horizontal intervals on the line y to leave open around circles and text."""

    out=[]

    for d in page.get_drawings():

        r=d['rect']

        circle=d['type']=='fs' and r.width<15

        halo=d.get('color') and tuple(round(c,3) for c in d['color'])==(1.0,1.0,1.0) and round(d.get('width') or 0,2)==1.7

        if (circle or halo) and r.y0-pad<=y<=r.y1+pad:out.append((r.x0-pad,r.x1+pad))

    for b in page.get_text('dict')['blocks']:

        for l in b.get('lines',[]):

            for s in l['spans']:

                x0,y0,x1,y1=s['bbox']

                if y0-pad<=y<=y1+pad and s['text'].strip():out.append((x0-pad,x1+pad))

    return sorted(out)



def segments(x0,x1,blocked):

    pieces=[];start=x0

    for a,b in blocked:

        if b<=start or a>=x1:continue

        if a>start:pieces.append((start,a))

        start=max(start,b)

    if start<x1:pieces.append((start,x1))

    return pieces



def main():

    doc=pymupdf.open(SRC);page=doc[0]

    rivers=strip_rivers(doc,page)

    frame=[d['rect'] for d in page.get_drawings() if d.get('color') and tuple(round(c,3) for c in d['color'])==(0.2,0.2,0.2)]

    fx0=min(r.x0 for r in frame);fx1=max(r.x1 for r in frame);fy0=min(r.y0 for r in frame);fy1=max(r.y1 for r in frame)

    ylat=lambda lat:fy1-(fy1-fy0)*lat/50.0

    xref=page.get_images()[0][0];info=doc.extract_image(xref)



    hires=output_path('data/figure_inputs/mff_relief_z7_platecarree.png')

    rgb=np.array(Image.open(hires).convert('RGB')) if hires.exists() else np.array(Image.open(io.BytesIO(info['image'])).convert('RGB'))

    new,sea_fraction=terrain_tint(rgb)

    buf=io.BytesIO();Image.fromarray(new).save(buf,format='JPEG',quality=92,subsampling=0)

    page.replace_image(xref,stream=buf.getvalue())

    # Text positions are read before drawing so the new lines leave room for existing labels.

    lines=[(40,(0,(4,2.5)),'Track-composition domain limit (40°N)'),(20,(0,(1,2)),'20°N')]

    for lat,_,_ in lines:pass

    for lat,dash,label in lines:

        y=ylat(lat);blocked=gaps(page,y)

        for a,b in segments(fx0,fx1,blocked):

            page.draw_line((a,y),(b,y),color=LINE_COLOR,width=0.9,dashes='[4 2.5] 0' if lat==40 else '[1 2] 0')

    page.insert_font(fontname='LibSans',fontfile=FONT)

    for lat,_,label in lines:

        y=ylat(lat);size=6.2

        w=pymupdf.Font(fontfile=FONT).text_length(label,fontsize=size)

        x=fx1-4-w

        page.draw_rect(pymupdf.Rect(x-1.5,y-size-2.2,x+w+1.5,y-1.2),color=None,fill=(1,1,1),fill_opacity=0.75)

        page.insert_text((x,y-2.6),label,fontname='LibSans',fontsize=size,color=LINE_COLOR)

    doc.set_metadata({**doc.metadata,'title':'Western North Pacific analysis domain and mutually exclusive coastal segments',

                      'subject':'Journal of Geophysical Research: Atmospheres Figure 1','creator':'Matplotlib','author':'Museman X'})

    # Write through Python: the mounted folder forbids the delete-then-write used by mupdf.

    pdf=output_path('figures/Fig01_study_region.pdf');data=doc.tobytes(garbage=3,deflate=True);pdf.write_bytes(data)

    pix=pymupdf.open(stream=data,filetype='pdf')[0].get_pixmap(dpi=600);output_path('figures/Fig01_study_region.png').write_bytes(pix.tobytes('png'))

    write_json('results/diagnostics/study_region.json',dict(status='PASS',basemap_pixels=list(new.shape[1::-1]),source='Reference PDF layout; Maps-For-Free relief at 30% colour saturation and 30% opacity over white',terrain_tint=TERRAIN_TINT,terrain_saturation=TERRAIN_SATURATION,lightening_method="PIL Color factor 0.30, then white blend at opacity 0.30",river_groups_removed=rivers,

        sea_pixel_fraction=round(float(sea_fraction),4),added=['40°N northern limit of track-composition domain (dashed)','20°N (dotted)'],png_size=[pix.width,pix.height]))

    print('Figure 1 written',pix.width,pix.height)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()
