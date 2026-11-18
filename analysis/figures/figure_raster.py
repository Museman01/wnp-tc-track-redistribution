"""Export figure PNGs and optionally resize an embedded terrain raster."""
from io import BytesIO
from pathlib import Path
from tempfile import NamedTemporaryFile

from PIL import Image
import pymupdf


def render_png(pdf, dpi=600):
    pdf = Path(pdf)
    with pymupdf.open(pdf) as document:
        page = document[0].get_pixmap(dpi=dpi, colorspace=pymupdf.csRGB, alpha=False)
        image = Image.frombytes('RGB', (page.width, page.height), page.samples)
    with pdf.with_suffix('.png').open('xb') as output:
        image.save(output, format='PNG', dpi=(dpi, dpi), optimize=dpi == 300)


def resize_terrain(pdf, size, scale, count, *, jpeg_quality=None):
    pdf = Path(pdf)
    with NamedTemporaryFile(prefix=pdf.stem + '-', suffix='.pdf', dir=pdf.parent, delete=False) as tmp:
        temporary = Path(tmp.name)
    try:
        with pymupdf.open(pdf) as document:
            targets = [item for item in document[0].get_images(full=True)
                       if all(abs(a-b) <= 1 for a,b in zip(item[2:4], size))]
            if len(targets) != count:
                raise ValueError(f'Expected {count} terrain rasters of size {size}: {pdf}')
            for xref, mask, width, height, *_ in targets:
                raw = document.xref_stream(xref)
                if mask or len(raw) != width * height * 3:
                    raise ValueError(f'Terrain raster must be RGB with no mask: {pdf}')
                image = Image.frombytes('RGB', (width, height), raw)
                image = image.resize(tuple(round(v * scale) for v in size), Image.Resampling.LANCZOS)
                document.xref_set_key(xref, 'Width', str(image.width))
                document.xref_set_key(xref, 'Height', str(image.height))
                document.xref_set_key(xref, 'DecodeParms', 'null')
                if jpeg_quality is None:
                    document.update_stream(xref, image.tobytes(), compress=True)
                else:
                    encoded = BytesIO()
                    image.save(encoded, 'JPEG', quality=jpeg_quality, subsampling=0, optimize=True)
                    document.update_stream(xref, encoded.getvalue(), compress=False)
                    document.xref_set_key(xref, 'Filter', '/DCTDecode')
            document.save(temporary, garbage=4, deflate=True)
        temporary.replace(pdf)
    finally:
        temporary.unlink(missing_ok=True)
