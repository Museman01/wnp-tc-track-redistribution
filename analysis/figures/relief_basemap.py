'Relief basemap.'
from __future__ import annotations

import sys,math,time,urllib.request

from concurrent.futures import ThreadPoolExecutor

import numpy as np

from PIL import Image

from core.paths import PROJECT,output_path



ZOOM=7

TILE_DIR=output_path(f'data/figure_inputs/mff_relief_z{ZOOM}')

OUT=output_path(f'data/figure_inputs/mff_relief_z{ZOOM}_platecarree.png')

URL='https://maps-for-free.com/layer/relief/z{z}/row{y}/{z}_{x}-{y}.jpg'

# Geographic extent of the raster as placed in the reference PDF (image box relative to the

# 100-160E, 0-50N axes frame).

EXTENT=dict(lon0=100.003,lon1=160.015,lat0=0.067,lat1=50.024)

WIDTH=5000



def tile_xy(lon,lat,z=ZOOM):

    n=2**z;x=(lon+180)/360*n;y=(1-math.asinh(math.tan(math.radians(lat)))/math.pi)/2*n

    return x,y



def tiles():

    x0,y0=tile_xy(EXTENT['lon0'],EXTENT['lat1']);x1,y1=tile_xy(EXTENT['lon1'],EXTENT['lat0'])

    return [(x,y) for x in range(int(x0),int(x1)+1) for y in range(int(y0),int(y1)+1)]



def fetch(xy):

    x,y=xy;path=TILE_DIR/f'{x}_{y}.jpg'

    if path.exists() and path.stat().st_size>0:return 'cached'

    for attempt in range(3):

        try:

            req=urllib.request.Request(URL.format(z=ZOOM,x=x,y=y),headers={'User-Agent':'research-figure/1.0'})

            data=urllib.request.urlopen(req,timeout=30).read()

            if data[:2]==b'\xff\xd8':path.write_bytes(data);return 'ok'

        except Exception:time.sleep(1+attempt)

    return 'failed'



def download(budget=150):

    TILE_DIR.mkdir(parents=True,exist_ok=True);todo=tiles();start=time.monotonic();status=[]

    with ThreadPoolExecutor(8) as pool:

        for i in range(0,len(todo),32):

            if time.monotonic()-start>budget:break

            status+=list(pool.map(fetch,todo[i:i+32]))

    have=sum((TILE_DIR/f'{x}_{y}.jpg').exists() for x,y in todo)

    print(f'tiles {have}/{len(todo)}',dict((s,status.count(s)) for s in set(status)))



def build(output):

    out=(PROJECT / output).resolve()
    if not out.is_relative_to(PROJECT.resolve()): raise ValueError('Output must remain within the project')
    if out.exists(): raise FileExistsError('Refusing to overwrite an existing basemap')
    out.parent.mkdir(parents=True,exist_ok=True)
    from scipy.ndimage import map_coordinates

    todo=tiles();xs=sorted({x for x,_ in todo});ys=sorted({y for _,y in todo})

    mosaic=np.zeros((len(ys)*256,len(xs)*256,3),np.uint8)

    for x,y in todo:

        mosaic[(y-ys[0])*256:(y-ys[0]+1)*256,(x-xs[0])*256:(x-xs[0]+1)*256]=np.array(Image.open(TILE_DIR/f'{x}_{y}.jpg').convert('RGB'))

    # Square pixels on the page: the PDF image box is 424.08 x 389.52 pt.

    height=round(WIDTH*389.52/424.08)

    lon=EXTENT['lon0']+(np.arange(WIDTH)+0.5)/WIDTH*(EXTENT['lon1']-EXTENT['lon0'])

    lat=EXTENT['lat1']-(np.arange(height)+0.5)/height*(EXTENT['lat1']-EXTENT['lat0'])

    px=(lon+180)/360*2**ZOOM*256-xs[0]*256

    py=(1-np.arcsinh(np.tan(np.radians(lat)))/np.pi)/2*2**ZOOM*256-ys[0]*256

    PY,PX=np.meshgrid(py,px,indexing='ij')

    out=np.stack([map_coordinates(mosaic[...,c].astype(float),[PY-0.5,PX-0.5],order=1,mode='nearest') for c in range(3)],-1)

    Image.fromarray(out.clip(0,255).astype(np.uint8)).save((PROJECT / output).resolve())

    print('basemap',WIDTH,height)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    {'download':download,'build':build}[sys.argv[1]](*sys.argv[2:])

