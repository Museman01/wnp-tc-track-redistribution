"""Center the existing seven map digits on their existing circular markers."""
from __future__ import annotations
from analysis.figures.figure_raster import render_png, resize_terrain
import argparse
import io
import re
from pathlib import Path

import pymupdf
from fontTools.ttLib import TTFont

PAT = re.compile(rb"q ([^qQ]*?)BT/(Retained\d+) ([\d.]+) Tf ([-\d.]+) ([-\d.]+) TD\[(.*?)\]TJ ET Q")


def align(source: Path, output: Path):
    if any(Path(str(output) + suffix).exists() for suffix in ('.pdf', '.png')):
        raise FileExistsError(output)
    doc = pymupdf.open(source)
    page = doc[0]
    circles = []
    for d in page.get_drawings():
        r = d['rect']
        if (d['type'] == 'fs' and d.get('fill') == (1.0, 1.0, 1.0)
                and 6 < r.width < 20 and 6 < r.height < 20 and abs(r.width-r.height) < 1):
            circles.append(r)
    if len(circles) != 7:
        raise ValueError(f'Expected seven circles, found {len(circles)}')
    bold_fonts = [x for x in doc.get_page_fonts(0, full=True) if x[3] == 'Arial Bold']
    if len(bold_fonts) != 1:
        raise ValueError('Expected the original embedded Arial Bold font')
    font = TTFont(io.BytesIO(doc.extract_font(bold_fonts[0][0])[3]))
    cmap, em, glyf = font.getBestCmap(), font['head'].unitsPerEm, font['glyf']
    chars = []
    for block in page.get_text('rawdict')['blocks']:
        for line in block.get('lines', []):
            for span in line['spans']:
                if len(span['chars']) != 1 or span['font'] != 'Arial-BoldMT':
                    continue
                c = span['chars'][0]
                if c['c'] not in '1234567':
                    continue
                x, y = c['origin']; size = span['size']
                circle = min(circles, key=lambda r: (x + size * 0.28 - (r.x0+r.x1)/2)**2 + (y-size*0.36-(r.y0+r.y1)/2)**2)
                glyph = glyf[cmap[ord(c['c'])]]
                ink_x = x + (glyph.xMin+glyph.xMax)*size/(2*em)
                ink_y = y - (glyph.yMin+glyph.yMax)*size/(2*em)
                dx = (circle.x0+circle.x1)/2 - ink_x
                dy = (circle.y0+circle.y1)/2 - ink_y
                chars.append((c['c'], x, y, dx, dy))
    if sorted(c[0] for c in chars) != list('1234567'):
        raise ValueError('Expected exactly one instance of each numbered marker')
    xref = page.get_contents()[0]
    data = doc.xref_stream(xref)
    changed = []
    def replace(match):
        prefix, fontname, size, lx, ly, body = match.groups()
        if fontname != bold_fonts[0][4].encode() or b' cm' in prefix:
            return match.group()
        x, y = float(lx), page.rect.height-float(ly)
        c = min(chars, key=lambda t: (t[1]-x)**2+(t[2]-y)**2)
        if (c[1]-x)**2+(c[2]-y)**2 > .01:
            return match.group()
        changed.append(c[0])
        nx, ny = x+c[3], y+c[4]
        return (b'q '+prefix+b'BT/'+fontname+b' '+size+f' Tf {nx:.6f} {page.rect.height-ny:.6f} TD['.encode()+body+b']TJ ET Q')
    data = PAT.sub(replace, data)
    if sorted(changed) != list('1234567'):
        raise ValueError(f'Unexpected glyph edits: {changed}')
    doc.update_stream(xref, data)
    output.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(output)+'.pdf', garbage=3, deflate=True)
    doc.close()
    pdf = Path(str(output)+'.pdf')
    resize_terrain(pdf, (5000, 4593), 0.54, 1, jpeg_quality=93)
    render_png(pdf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    align(args.input, args.output)
