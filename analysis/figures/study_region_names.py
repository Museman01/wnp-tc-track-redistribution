from pathlib import Path
import re
import pymupdf

def strip(source,out):
 if out.exists() or (out.parent/"sea_name_bounds.txt").exists():raise FileExistsError(out)
 doc=pymupdf.open(source);page=doc[0];drawings=[x for x in page.get_drawings() if x.get("fill") and abs(x["fill"][0]-.08235294)<1e-4 and abs(x["fill"][1]-.352941)<1e-4]
 if len(drawings)!=20:raise ValueError(len(drawings))
 words=[drawings[i]["rect"] for i in range(0,20,2)]
 (out.parent/"sea_name_bounds.txt").write_text("\n".join(" ".join(str(v) for v in r) for r in words),encoding="utf-8")
 xref=page.get_contents()[0];data=doc.xref_stream(xref);start=data.find(b"/A9 gs .08235294 .3529412 .57254907 rg")
 if start<0:raise ValueError("marker")
 pat=re.compile(rb"[-.0-9]+ [-.0-9]+ m ");endpat=re.compile(rb" (?:B|f) ");pos=start;edits=[]
 for index in range(20):
  a=pat.search(data,pos);b=endpat.search(data,a.end()) if a else None
  if not b:raise ValueError(index)
  edits.append((a.start(),b.end()));pos=b.end()
 for a,b in reversed(edits):data=data[:a]+data[b:]
 doc.update_stream(xref,data);doc.save(out,garbage=3,deflate=True);doc.close()
def draw_names(source,bounds,out):
 from matplotlib.textpath import TextPath
 from matplotlib.font_manager import FontProperties
 from matplotlib.path import Path as MPath
 from analysis.figures.figure_raster import render_png
 doc=pymupdf.open(source);page=doc[0];path="C:/Windows/Fonts/timesi.ttf";font=pymupdf.Font(fontfile=path);prop=FontProperties(fname=path);page.insert_font(fontname="SeaTNR",fontfile=path)
 names=["South China","Sea","East China","Sea","Sea of","Japan","Philippine","Sea","Western North","Pacific"]
 for name,line in zip(names,Path(bounds).read_text(encoding="utf-8").splitlines()):
  rect=pymupdf.Rect([float(x) for x in line.split()]);rect.y0+=18;rect.y1+=18;size=12.4;w=font.text_length(name,fontsize=size);x=(rect.x0+rect.x1-w)/2;y=rect.y1-.8
  glyph=TextPath((0,0),name,size=size,prop=prop);shape=page.new_shape();current=None;start=None
  def pt(a,b):return pymupdf.Point(x+float(a),y-float(b))
  for coords,code in glyph.iter_segments(curves=True):
   if code==MPath.MOVETO:current=pt(*coords);start=current
   elif code==MPath.LINETO:
    q=pt(*coords);shape.draw_line(current,q);current=q
   elif code==MPath.CURVE3:
    ctrl=pt(*coords[:2]);end=pt(*coords[2:]);c1=current+(ctrl-current)*(2/3);c2=end+(ctrl-end)*(2/3);shape.draw_bezier(current,c1,c2,end);current=end
   elif code==MPath.CURVE4:
    c1=pt(*coords[:2]);c2=pt(*coords[2:4]);end=pt(*coords[4:]);shape.draw_bezier(current,c1,c2,end);current=end
   elif code==MPath.CLOSEPOLY and current is not None and start is not None:shape.draw_line(current,start);current=start
  shape.finish(color=(1,1,1),width=.8,closePath=False);shape.commit();page.insert_text((x,y),name,fontname="SeaTNR",fontsize=size,color=(.08235,.35294,.57255))
 out=Path(out)
 if out.exists() or out.with_suffix(".png").exists():raise FileExistsError(out)
 out.parent.mkdir(parents=True,exist_ok=True);doc.save(out,garbage=3,deflate=True);doc.close();render_png(out)
if __name__=="__main__":
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument("mode",choices=["strip","names"]);parser.add_argument("--input",required=True);parser.add_argument("--output",required=True,type=Path);parser.add_argument("--bounds",type=Path);args=parser.parse_args()
 if args.mode=="strip":args.output.parent.mkdir(parents=True,exist_ok=True);strip(args.input,args.output)
 else:
  if not args.bounds:parser.error("names requires --bounds")
  draw_names(args.input,args.bounds,args.output)
