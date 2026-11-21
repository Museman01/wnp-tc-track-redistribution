from pathlib import Path
import re
import pymupdf
from analysis.figures.retained_text import PAT
from analysis.figures.figure_raster import render_png

def main(source, output):
    source=Path(source)
    out=Path(str(output)+".pdf")
    if out.exists() or out.with_suffix(".png").exists(): raise FileExistsError(out)
    doc=pymupdf.open(source);page=doc[0]; selected=[]
    for block in page.get_text("rawdict")["blocks"]:
     for line in block.get("lines",[]):
      chars=[c for span in line["spans"] for c in span["chars"]]
      text="".join(c["c"] for c in chars)
      is_longitude=line["bbox"][1]>400 and text.endswith("E")
      is_scale=line["bbox"][0]>340 and line["bbox"][1]>350
      if is_longitude or is_scale:
       selected.extend((c["origin"][0],c["origin"][1],c["c"]) for c in chars)
    xref=page.get_contents()[0];data=doc.xref_stream(xref); removed=[]
    def strip(m):
     prefix,font,size,lx,ly,body=m.groups();x,y=float(lx),page.rect.height-float(ly)
     if b" cm" in prefix:return m.group()
     if not selected:return m.group()
     near=min(selected,key=lambda v:(v[0]-x)**2+(v[1]-y)**2)
     if (near[0]-x)**2+(near[1]-y)**2<.01:
      removed.append(near[2]);return b""
     return m.group()
    data=PAT.sub(strip,data)
    assert len(removed)==len(selected),(len(removed),len(selected))
    for x in ("33.81875","175.14987","316.481","457.81214"):
     pat=rf"{x} 20\.360156 m {x} 16\.860156 l B".encode();data,n=re.subn(pat,b"",data);assert n==1,(x,n)
    data=data.replace(b"349.34049 42.191987 m 446.0289 42.191987 l",b"349.34049 42.191987 m 413.79943 42.191987 l")
    data=data.replace(b"397.6847 39.540979 m 397.6847 44.842996 l",b"381.56996 39.540979 m 381.56996 44.842996 l")
    data=data.replace(b"446.0289 39.540979 m 446.0289 44.842996 l",b"413.79943 39.540979 m 413.79943 44.842996 l")
    doc.update_stream(xref,data)
    output=pymupdf.open();pg=output.new_page(width=page.rect.width,height=419)
    pg.show_pdf_page(pymupdf.Rect(0,18,page.rect.width,418),doc,0,clip=pymupdf.Rect(0,0,page.rect.width,400))
    fontpath="C:/Windows/Fonts/arial.ttf";pg.insert_font(fontname="ArialUser",fontfile=fontpath);font=pymupdf.Font(fontfile=fontpath)
    def text_center(x,y,text,size=10.2):
     pg.insert_text((x-font.text_length(text,fontsize=size)/2,y),text,fontname="ArialUser",fontsize=size,color=(.125,.129,.141))
    for x,label in zip([33.81875,175.14987,316.48099,457.81214],["100°E","120°E","140°E","160°E"]):
     pg.draw_line((x,25.385),(x,21.885),color=(.15,.15,.15),width=.8);text_center(x,18,label)
    text_center(349.34049,407.0,"0");text_center(381.56996,388.8,"500");text_center(413.79943,407.0,"1000 km")
    out.parent.mkdir(parents=True,exist_ok=True);output.save(out,garbage=3,deflate=True);output.close();doc.close();render_png(out)
    print("Removed glyphs",len(removed),"and wrote",out)



if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    main(args.input,args.output)
