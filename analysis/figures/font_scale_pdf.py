from analysis.figures.figure_raster import render_png
from pathlib import Path
import pymupdf
from analysis.figures import retained_text
import argparse

def run(source,output):
 d=pymupdf.open(source);retained_text.enlarge(d,factor=1.42,extra_multiline_gap=6)
 output=Path(output)
 if any(Path(str(output)+e).exists() for e in ['.pdf','.png']):raise FileExistsError(output)
 output.parent.mkdir(parents=True,exist_ok=True);d.save(str(output)+'.pdf',garbage=3,deflate=True)
 render_png(str(output)+'.pdf')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args();run(a.input,a.output)
