"""Enlarge individually positioned PDF text while scaling its character spacing."""
from __future__ import annotations
import re

PAT = re.compile(rb"q ([^qQ]*?)BT/(Retained\d+) ([\d.]+) Tf ([-\d.]+) ([-\d.]+) TD\[(.*?)\]TJ ET Q")

def enlarge(doc, factor=1.10, extra_multiline_gap=0, preserve_marker_digits=False, right_axis_shift=0):
    page=doc[0]; chars=[]
    for block in page.get_text('rawdict')['blocks']:
        for line in block.get('lines', []):
            entries=[c for span in line['spans'] for c in span['chars']]
            if not entries: continue
            dx,dy=line['dir']; projected=[c['origin'][0]*dx+c['origin'][1]*dy for c in entries]
            center=(min(projected)+max(projected))/2
            text=''.join(c['c'] for c in entries)
            gap=extra_multiline_gap if text.startswith('vs') else 0
            shift=right_axis_shift if abs(dx)<0.01 and line['bbox'][0]>page.rect.width*0.8 and len(text)>4 else 0
            local_factor=1.0 if preserve_marker_digits and text in '1234567' and len(text)==1 else factor
            for c,projection in zip(entries,projected):
                x,y=c['origin']; delta=(local_factor-1)*(projection-center)
                chars.append((x,y,x+delta*dx+shift,y+delta*dy+gap,local_factor))
    xref=page.get_contents()[0]; data=doc.xref_stream(xref); found=0
    def change(m):
        nonlocal found
        prefix,font,size,lx,ly,body=m.groups(); x,y=float(lx),float(ly); rotated=b'0 1 -1 0 0 0 cm' in prefix
        px,py=(-y,page.rect.height-x) if rotated else (x,page.rect.height-y)
        near=min(chars,key=lambda t:(t[0]-px)**2+(t[1]-py)**2)
        if (near[0]-px)**2+(near[1]-py)**2 > .01: raise ValueError('Unmapped PDF glyph')
        nx,ny=near[2],near[3]
        ax,ay=(page.rect.height-ny,-nx) if rotated else (nx,page.rect.height-ny)
        size=float(size)*near[4]; found+=1
        return b'q '+prefix+b'BT/'+font+f' {size:.5f} Tf {ax:.6f} {ay:.6f} TD['.encode()+body+b']TJ ET Q'
    data=PAT.sub(change,data)
    if found!=len(chars):raise ValueError(f'PDF glyph count differs: {found} vs {len(chars)}')
    doc.update_stream(xref,data)
    return found
