"""Enlarge PDF text while retaining the existing scientific vector and raster layers."""
from pathlib import Path
import json, hashlib, re
from collections import Counter
import pymupdf as fitz
import numpy as np

from core.paths import PROJECT_ROOT as PROJECT
INPUT=PROJECT/'data/figure_inputs/typography'
PARAMETERS=json.loads((INPUT/'display_parameters.json').read_text())
OUTPUT=None
CONFIG={'Fig01_study_region':9.,'FigS04_activity_timeseries':11.5,'FigS05_lmi_location':11.,'FigS08_wnpsh_fixed_contour_metrics':11.5,'FigS10_ohc_latitude_contributions':11.5}

def normalized_drawings(page):
    return [{k:v for k,v in x.items() if k not in ['seqno','layer']} for x in page.get_drawings()
            if not (len(x['items'])==2 and 70<x['rect'].x0<84 and x['rect'].x1<97 and x['rect'].height<.01 and 23<x['rect'].y0<50)
            and not (x['type']=='f' and 373<x['rect'].x0<386 and 11<x['rect'].y0<12 and x['rect'].width<10 and x['rect'].height<7)
            and not (page.parent.metadata.get('title','').startswith('Western North Pacific analysis domain') and x.get('fill')==(1.,1.,1.) and abs(x.get('fill_opacity',1)-.75)<.01)]

def redraw(stem):
    doc=fitz.open(INPUT/(stem+'.pdf'));page=doc[0]
    before_draw=normalized_drawings(page)
    before_images=[hashlib.sha256(doc.extract_image(x[0])['image']).hexdigest() for x in page.get_images()]
    lines=[line for block in page.get_text('rawdict')['blocks'] for line in block.get('lines',[])]
    fonts={}
    for n,row in enumerate(page.get_fonts()):
        data=doc.extract_font(row[0])[3]
        base_name=row[3].split('+')[-1]
        if base_name.startswith('Arial'):
            filename='arialbi.ttf' if 'BoldItalic' in base_name else 'arialbd.ttf' if 'Bold' in base_name else 'ariali.ttf' if 'Italic' in base_name else 'arial.ttf'
            font_path=Path('C:/Windows/Fonts')/filename
            if font_path.is_file():data=font_path.read_bytes()
        if not data:raise ValueError(row)
        base=row[3].split('+')[-1];tag='Retained'+str(n)
        page.insert_font(fontname=tag,fontbuffer=data)
        fonts[base]=(tag,fitz.Font(fontbuffer=data),data)
    original_characters=Counter(c['c'].replace('\xad','-') for line in lines for span in line['spans'] for c in span['chars'] if not c['c'].isspace())
    records=[];characters=[];backgrounds=[]
    for line in lines:
        text=''.join(c['c'] for span in line['spans'] for c in span['chars'])
        body=max(s['size'] for s in line['spans']);target=CONFIG[stem]
        left_anchored=text in ['Annual value','Sen trend','Trend','1966–1981','1982–2025','South of 20°N','20°N and north']
        # One ordinary size for axes, legends, annotations and panel letters.
        if stem=='Fig01_study_region':target=10.2
        scale=target/body
        direction=np.array(line['dir']);start=np.array(line['spans'][0]['origin'])
        bbox=fitz.Rect(line['bbox'])
        length=bbox.width if abs(direction[0])>.5 else bbox.height
        anchor=start if left_anchored else start+direction*length/2
        if stem=='Fig01_study_region' and ('domain' in text or body<7):anchor=start+direction*length
        rotate=90 if direction[1]<-.5 else 270 if direction[1]>.5 else 180 if direction[0]<-.5 else 0
        record=dict(text=text,target=target,source_size=body,rotate=rotate)
        records.append(record)
        if stem=='Fig01_study_region' and 'domain' in text:
            # Reflow the same annotation into open water beneath its 40-degree line.
            span=line['spans'][0];tag,font,data=fonts[span['font']];right=bbox.x1
            for row,label in enumerate(['Track-composition','domain limit','(40°N)']):
                width=font.text_length(label,fontsize=target);x=right-width;y=start[1]+target*(2.3+row*1.15)
                backgrounds.append(fitz.Rect(x-1.5,y-target-1.5,right+1.5,y+2))
                for ch in label:
                    characters.append((ch,np.array([x,y]),tag,target,fitz.sRGB_to_pdf(span['color']),0));x+=font.text_length(ch,fontsize=target)
            continue
        for span in line['spans']:
            fontname=span['font'];tag,font,data=fonts[fontname]
            size=span['size']*scale
            color=fitz.sRGB_to_pdf(span['color'])
            for char in span['chars']:
                value=char['c'];old=np.array(char['origin']);origin=anchor+(old-anchor)*scale
                if stem=='FigS08_wnpsh_fixed_contour_metrics' and left_anchored:origin[0]-=13
                if stem=='FigS10_ohc_latitude_contributions' and text=='20°N and north':origin[0]+=12
                characters.append((value,origin,tag,size,color,rotate))
    # Remove only PDF text objects; preserve every non-text content-stream byte.
    stream_checks=[]
    for xref in page.get_contents():
        original=doc.xref_stream(xref)
        if stem=='Fig01_study_region' and len(original)<180 and b'1 1 1 rg f' in original and b' gs' in original:
            doc.update_stream(xref,b' ');continue
        if stem=='FigS08_wnpsh_fixed_contour_metrics':
            def shift_handle(m):return f"{float(m[1])-13:.6f} {m[2].decode()} {m[3].decode()}".encode()
            original,count_handles=re.subn(rb'(?m)^(83\.424474|89\.864474|96\.304474) ([0-9.]+) ([ml])$',shift_handle,original)
            assert count_handles in (0,9),('legend handles',count_handles)
        if stem=='FigS10_ohc_latitude_contributions':
            def shift_second_handle(m):return f"{float(m[1])+12:.6f} {m[2].decode()} {m[3].decode()}".encode()
            original,n=re.subn(rb'(?m)^(373\.235073|382\.435072) (273\.765908|280\.205908) ([ml])$',shift_second_handle,original)
            assert n in (0,4),n
        changed,count=re.subn(rb'(?s)(?<!\S)BT(?!\S).*?(?<!\S)ET(?!\S)',b'',original)
        if count:doc.update_stream(xref,changed)
        stream_checks.append([xref,count])
    assert not page.get_text().strip(),(stem,'text not fully removed')
    assert normalized_drawings(page)==before_draw,(stem,'vector layer changed')
    for tag,font,data in fonts.values():page.insert_font(fontname=tag,fontbuffer=data)
    for rect in backgrounds:page.draw_rect(rect,fill=(1,1,1),color=None,fill_opacity=.75)
    for text,origin,tag,size,color,rotate in characters:
        page.insert_text(fitz.Point(*origin),text,fontname=tag,fontsize=size,color=color,rotate=rotate,overlay=True)
    assert normalized_drawings(page)==before_draw
    assert [hashlib.sha256(doc.extract_image(x[0])['image']).hexdigest() for x in page.get_images()]==before_images
    # Reinserted mathematical superscripts keep their relative scale and position.
    expected=Counter(c[0].replace('\xad','-') for c in characters if not c[0].isspace())
    assert expected==original_characters,(stem,'original characters changed')
    actual=Counter(c.replace('\xad','-') for c in page.get_text() if not c.isspace())
    assert expected==actual,(stem,expected-actual,actual-expected)
    output=OUTPUT/(stem+'.pdf');doc.save(output,garbage=3,deflate=True)
    from PIL import Image
    export_dpi=PARAMETERS[stem]['dpi']
    page.get_pixmap(dpi=export_dpi).save(OUTPUT/(stem+'.png'))
    print(stem,flush=True)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('figure', choices=sorted(CONFIG))
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if not args.execute:parser.error('Use --execute to draw figures')
    OUTPUT=args.output_dir.resolve();OUTPUT.mkdir(parents=True,exist_ok=True)
    if any(OUTPUT.iterdir()):parser.error('Select an empty output directory')
    redraw(args.figure)
