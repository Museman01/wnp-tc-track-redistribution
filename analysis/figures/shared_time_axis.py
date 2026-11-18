from analysis.figures.figure_raster import render_png
from pathlib import Path
import re, argparse
import pymupdf
from matplotlib.font_manager import FontProperties, findfont
from analysis.figures import retained_text

def render(source, output):
    output=Path(output)
    if any(Path(str(output)+ext).exists() for ext in ('.png','.pdf')): raise FileExistsError(output)
    src=pymupdf.open(source); p=src[0]; page_height=p.rect.height
    lines=[sp for b in p.get_text('dict')['blocks'] for ln in b.get('lines',[]) for sp in ln['spans']]
    years=[sp for sp in lines if re.fullmatch(r'(?:19|20)\d\d',sp['text']) and sp['bbox'][1] > page_height-45]
    if len(years) not in (12,14): raise ValueError(f'Unexpected bottom year labels: {len(years)}')
    cut=min(sp['bbox'][1] for sp in years)-1.0
    retained_text.enlarge(src,factor=1.39,right_axis_shift=9)
    bottom_baselines=[sp['origin'][1] for sp in years]
    stream=p.get_contents()[0]
    def strip_old_years(match):
        prefix, _, _, _, y, _ = match.groups()
        is_bottom_year = b'0 1 -1 0 0 0 cm' not in prefix and any(
            abs(page_height-float(y)-baseline)<0.05 for baseline in bottom_baselines)
        return b'' if is_bottom_year else match.group()
    src.update_stream(stream,retained_text.PAT.sub(strip_old_years,src.xref_stream(stream)))
    margin_x, margin_top = 12, 43
    dst=pymupdf.open(); page=dst.new_page(width=p.rect.width+2*margin_x,height=cut+margin_top+8)
    page.show_pdf_page(pymupdf.Rect(margin_x,margin_top,margin_x+p.rect.width,margin_top+cut),src,0,clip=pymupdf.Rect(0,0,p.rect.width,cut))
    if len(years) == 12:
        left = pymupdf.Rect(0,cut,49.0,min(cut+6,p.rect.height))
        page.show_pdf_page(pymupdf.Rect(margin_x,margin_top+cut,margin_x+left.width,margin_top+cut+left.height),src,0,clip=left)
    fontfile = Path(findfont(FontProperties(family=['Arial', 'Liberation Sans']), fallback_to_default=False)).as_posix()
    page.insert_font(fontname='LabelArial', fontfile=fontfile)
    font = pymupdf.Font(fontfile=fontfile)
    for sp in years:
        cx=margin_x+(sp['bbox'][0]+sp['bbox'][2])/2
        text=sp['text'];w=font.text_length(text,fontsize=16.0)
        if not (len(years) == 14 and text == '1960'):
            page.insert_text((cx-w/2,margin_top+1),text,fontname='LabelArial',fontsize=16.0)
        page.draw_line((cx,margin_top+4),(cx,margin_top+7),color=(0,0,0),width=.8)
    for group in (years[:len(years)//2],years[len(years)//2:]):
        cx=margin_x+(group[0]['bbox'][0]+group[-1]['bbox'][2])/2
        w=font.text_length('Year',fontsize=16.0)
        page.insert_text((cx-w/2,16),'Year',fontname='LabelArial',fontsize=16.0)
    output.parent.mkdir(parents=True,exist_ok=True);dst.save(str(output)+'.pdf',garbage=3,deflate=True)

    render_png(str(output)+'.pdf')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True);parser.add_argument('--output',required=True);a=parser.parse_args();render(a.input,a.output)
