"""Render Figure S16 with the caption-defined steering-vector scale."""
from pathlib import Path
from core.paths import PROJECT
from analysis.figures.diagnostic_figures import steering


def render(output_dir):
    out=Path(output_dir).resolve()
    if not out.is_relative_to(PROJECT.resolve()):
        raise ValueError('Output directory must remain within the project')
    if any((out/('FigS16_steering_layer_regressions'+ext)).exists() for ext in ('.png','.pdf')):
        raise FileExistsError(out)
    out.mkdir(parents=True,exist_ok=True)
    steering(output_dir=out,reference_vector_degrees=9.0)


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir',type=Path)
    render(parser.parse_args().output_dir)


if __name__=='__main__':main()
