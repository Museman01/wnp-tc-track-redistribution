from analysis.figures.figure_raster import render_png
from pathlib import Path
import re, argparse
import pymupdf
from analysis.figures import retained_text

def render(source, output):
    output=Path(output)
    if any(Path(str(output)+ext).exists() for ext in ('.png','.pdf')): raise FileExistsError(output)
    src=pymupdf.open(source); p=src[0]; page_height=p.rect.height
    lines=[sp for b in p.get_text('dict')['blocks'] for ln in b.get('lines',[]) for sp in ln['spans']]
    years=[sp for sp in lines if re.fullmatch(r'(?:19|20)\d\d',sp['text']) and sp['bbox'][1] > page_height-45]
    if len(years) not in (12,14): raise ValueError(f'Unexpected bottom year labels: {len(years)}')
    cut=min(sp['bbox'][1] for sp in years)-1.0
    # Simplify left-axis labels and move b/d panel letters to the right.
    removable = {'17.5', '22.5', '27.5', '32.5', '70', '80', '90', '100'}
    remove_origins, move_origins, numeric_origins, axis_labels, right_titles = [], [], [], [], []
    for block in p.get_text('rawdict')['blocks']:
        for line in block.get('lines', []):
            chars = [ch for span in line['spans'] for ch in span['chars']]
            content = ''.join(ch['c'] for ch in chars)
            if chars and content in removable and line['bbox'][0] < 50:
                remove_origins.extend(ch['origin'] for ch in chars)
            if chars and content in ('b', 'd') and 340 < line['bbox'][0] < 365:
                move_origins.extend(ch['origin'] for ch in chars)
            if chars and abs(line['dir'][0]) < .01 and line['bbox'][0] > 650:
                right_titles.extend(ch['origin'] for ch in chars)
            if chars and re.fullmatch(r'\d+(?:\.\d+)?', content):
                x0, y0, x1, y1 = line['bbox']
                if y0 < page_height - 45 and (x0 < 50 or x0 > 640):
                    numeric_origins.extend(ch['origin'] for ch in chars)
                    if content not in removable or x0 > 640:
                        axis_labels.append((str(float(content)).removesuffix('.0'), x0 < 50, chars[0]['origin'][1], x0))
    assert len(remove_origins) == 25 and len(move_origins) == 2
    xref = p.get_contents()[0]
    def adjust_glyph(match):
        prefix, font, size, x, y, body = match.groups()
        rotated = b'0 1 -1 0 0 0 cm' in prefix
        pos = (-float(y), page_height-float(x)) if rotated else (float(x), page_height-float(y))
        near = lambda origin: (origin[0] - pos[0])**2 + (origin[1] - pos[1])**2 < .01
        if rotated and any(near(origin) for origin in right_titles):
            return b''
        if any(near(origin) for origin in remove_origins):
            return b''
        if any(near(origin) for origin in move_origins):
            return b''
        if any(near(origin) for origin in numeric_origins):
            return b''
        return match.group()
    src.update_stream(xref, retained_text.PAT.sub(adjust_glyph, src.xref_stream(xref)))
    retained_text.enlarge(src,factor=1.39,right_axis_shift=9)

    bottom_baselines=[sp['origin'][1] for sp in years]
    stream=p.get_contents()[0]
    def strip_old_years(match):
        prefix, _, _, _, y, _ = match.groups()
        is_bottom_year = b'0 1 -1 0 0 0 cm' not in prefix and any(
            abs(page_height-float(y)-baseline)<0.05 for baseline in bottom_baselines)
        return b'' if is_bottom_year else match.group()
    src.update_stream(stream,retained_text.PAT.sub(strip_old_years,src.xref_stream(stream)))
    margin_x, margin_top = 20, 43
    dst=pymupdf.open(); page=dst.new_page(width=p.rect.width+2*margin_x,height=cut+margin_top+8)
    page.show_pdf_page(pymupdf.Rect(margin_x,margin_top,margin_x+p.rect.width,margin_top+cut),src,0,clip=pymupdf.Rect(0,0,p.rect.width,cut))
    if len(years) == 12:
        left = pymupdf.Rect(0,cut,49.0,min(cut+6,p.rect.height))
        page.show_pdf_page(pymupdf.Rect(margin_x,margin_top+cut,margin_x+left.width,margin_top+cut+left.height),src,0,clip=left)
    page.insert_font(fontname='PanelBold',fontfile='C:/Windows/Fonts/arialbd.ttf')
    for label, y in (('b', 18.3), ('d', 208.5)):
        page.insert_text((margin_x + 621, margin_top + y), label, fontname='PanelBold', fontsize=16)
    fontfile='C:/Windows/Fonts/arial.ttf'; page.insert_font(fontname='LabelArial',fontfile=fontfile); font=pymupdf.Font(fontfile=fontfile)
    for label, middle in (('Storms yr', 101.5), ('ACE units yr', 287.0)):
        size = 16.0
        exponent_size = 11.0
        total = font.text_length(label, fontsize=size) + font.text_length('−1', fontsize=exponent_size)
        y = margin_top + middle + total / 2
        page.insert_text((margin_x + 690, y), label, fontname='LabelArial', fontsize=size, rotate=90)
        page.insert_text((margin_x + 686, y - font.text_length(label, fontsize=size)),
                         '−1', fontname='LabelArial', fontsize=exponent_size, rotate=90)
    for text, is_left, baseline, original_x in axis_labels:
        size = 16.0
        x = margin_x + (47.0 - font.text_length(text, fontsize=size) if is_left else original_x + 1.0)
        page.insert_text((x, margin_top + baseline), text, fontname='LabelArial', fontsize=size)
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
    print(source, 'copied years:',len(years),'cut:',cut)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True);parser.add_argument('--output',required=True);a=parser.parse_args();render(a.input,a.output)

