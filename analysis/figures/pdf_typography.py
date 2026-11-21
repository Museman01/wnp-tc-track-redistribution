from pathlib import Path
import argparse
import pymupdf
from analysis.figures.figure_raster import render_png

def convert(src,out):
 out=Path(out)
 if any(Path(str(out)+x).exists() for x in (".pdf",".png")):raise FileExistsError(out)
 doc=pymupdf.open(src);page=doc[0];spans=[];lines=[]
 for block in page.get_text("dict")["blocks"]:
  for line in block.get("lines",[]):
   group=[span for span in line["spans"] if span["text"].strip()]
   if group:lines.append((group,line["dir"],pymupdf.Rect(line["bbox"])))
   for span in group:spans.append((span,line["dir"]))
 for span,direction in spans:page.add_redact_annot(pymupdf.Rect(span["bbox"]),fill=False,cross_out=False)
 page.apply_redactions(images=0,graphics=0,text=0)
 fontpaths={"tnr":"C:/Windows/Fonts/times.ttf","tnrb":"C:/Windows/Fonts/timesbd.ttf","tnri":"C:/Windows/Fonts/timesi.ttf","tnrbi":"C:/Windows/Fonts/timesbi.ttf"}
 fonts={}
 for name,path in fontpaths.items():page.insert_font(fontname=name,fontfile=path);fonts[name]=pymupdf.Font(fontfile=path)
 for group,direction,bounds in lines:
  pieces=[]
  for span in group:
   name="tnr"+("b" if "Bold" in span["font"] else "")+("i" if "Italic" in span["font"] or "Oblique" in span["font"] else "")
   text=span["text"].replace("\xa0"," ");size=span["size"]*.9;w=fonts[name].text_length(text,fontsize=size);pieces.append((span,name,text,size,w))
  gaps=[]
  for left,right in zip(pieces,pieces[1:]):
   a,b=pymupdf.Rect(left[0]["bbox"]),pymupdf.Rect(right[0]["bbox"])
   gaps.append(max(0,(b.x0-a.x1) if direction[0]>.99 else (a.y0-b.y1))*.9)
  total=sum(piece[4] for piece in pieces)+sum(gaps)
  cursor=(bounds.x0+bounds.x1-total)/2 if direction[0]>.99 else (bounds.y0+bounds.y1+total)/2
  for i,(span,name,text,size,w) in enumerate(pieces):
   x,y=span["origin"];color=pymupdf.sRGB_to_pdf(span["color"])
   if direction[0]>.99:x=cursor;rotate=0;cursor+=w+(gaps[i] if i<len(gaps) else 0)
   elif direction[1]<-.99:y=cursor;rotate=90;cursor-=w+(gaps[i] if i<len(gaps) else 0)
   else:raise ValueError((text,direction))
   page.insert_text((x,y),text,fontname=name,fontsize=size,rotate=rotate,color=color)
 out.parent.mkdir(parents=True,exist_ok=True);doc.save(str(out)+".pdf",garbage=3,deflate=True);doc.close();render_png(str(out)+".pdf")
 print("TNR spans",len(spans))
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();convert(a.input,a.output)
