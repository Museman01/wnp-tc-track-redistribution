from pathlib import Path
import re
import pymupdf
from analysis.figures import retained_text

SCALE = 1.25

def main(source, output):
    SOURCE, OUT = Path(source), Path(output)
    if any(Path(str(OUT)+ext).exists() for ext in ('.pdf', '.png')):
        raise FileExistsError(OUT)
    doc = pymupdf.open(SOURCE); page = doc[0]
    drawings = [x for x in page.get_drawings() if x.get('fill') and abs(x['fill'][0] - 0.08235294) < 1e-4 and abs(x['fill'][1] - 0.352941) < 1e-4]
    if len(drawings) != 20:
        raise ValueError(f'Expected 20 sea-name outline paths, found {len(drawings)}')
    words = [drawings[i]['rect'] for i in range(0, 20, 2)]
    centers=[]
    for i in range(0,10,2):
        bounds=words[i] | words[i+1]
        centers.extend([(0.5*(bounds.x0+bounds.x1), page.rect.height-0.5*(bounds.y0+bounds.y1))]*4)
    retained_text.enlarge(doc,preserve_marker_digits=True)
    xref=page.get_contents()[0]; data=doc.xref_stream(xref)
    start=data.find(b'/A9 gs .08235294 .3529412 .57254907 rg')
    if start < 0: raise ValueError('Sea-name stream marker not found')
    path_start=re.compile(rb'[-.0-9]+ [-.0-9]+ m '); end_pattern=re.compile(rb' (?:B|f) ')
    pos=start; edits=[]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    for index,(cx,cy) in enumerate(centers):
        match=path_start.search(data,pos)
        end=end_pattern.search(data,match.end()) if match else None
        if not end: raise ValueError('Incomplete outline path')
        x,y=map(float,match.group()[:-3].split())
        if not (0 <= x <= page.rect.width and 0 <= y <= page.rect.height): raise ValueError('Unexpected path')
        matrix=f'q {SCALE} 0 0 {SCALE} {(1-SCALE)*cx - (5 if 16 <= index < 20 else 0):.6f} {(1-SCALE)*cy + (9 if 8 <= index < 12 else 0):.6f} cm '.encode()
        edits.append((match.start(),end.end(),matrix+data[match.start():end.end()]+b'Q '))
        pos=end.end()
    for a,b,v in reversed(edits): data=data[:a]+v+data[b:]
    doc.update_stream(xref,data)
    doc.save(str(OUT)+'.pdf',garbage=3,deflate=True)
    rendered=pymupdf.open(str(OUT)+'.pdf');pix=rendered[0].get_pixmap(dpi=600);pix.save(str(OUT)+'.png')
    print('Enlarged five ocean labels and the existing text; plot geometry was preserved.')
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    main(args.input,args.output)
