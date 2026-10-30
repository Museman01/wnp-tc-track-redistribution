"""Scientific analysis entry points."""
import argparse
import runpy
import sys

ROUTES = {
    "tracks": ("analysis.tracks.agency_comparison", []),
    "landfall": ("analysis.tracks.agency_comparison", ["landfall"]),
    "genesis": ("analysis.tracks.genesis_conditioned_decomposition", []),
    "speed": ("analysis.tracks.agency_comparison", ["speed"]),
    "steering": ("analysis.circulation.steering_layer_diagnostics", []),
    "trajectories": ("analysis.models.bam_primary_fit", []),
    "ocean": ("analysis.ocean.ohc_pipeline", []),
    "d26": ("analysis.ocean.d26_exposure", []),
    "pathways": ("analysis.tracks.pathway_closure", []),
    "figures": ("analysis.figures.main_figures", []),
    "diagnostics": ("analysis.figures.diagnostic_figures", []),
}

def main():
    parser = argparse.ArgumentParser(description="Track, landfall, circulation and ocean analyses.")
    parser.add_argument("analysis", choices=ROUTES)
    args, remaining = parser.parse_known_args()
    module, prefix = ROUTES[args.analysis]
    sys.argv = [module, *prefix, *remaining]
    runpy.run_module(module, run_name="__main__")

if __name__ == "__main__":
    main()
