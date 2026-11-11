"""Scientific analysis entry points."""
import argparse
import runpy
import sys

ROUTES = {
    'hko-robustness': ('analysis.tracks.hko_robustness_matrix', []),
    'regime-association': ('analysis.circulation.regime_recurvature_association', []),
    'regime-persistence': ('analysis.circulation.regime_persistence_recurvature', []),
    'height-gradient': ('analysis.circulation.height_gradient_diagnostic', []),
    'regime-logistic': ('analysis.circulation.regime_recurvature', []),
    'six-hourly-trajectories': ('analysis.models.six_hourly_trajectory', []),
    'submonthly': ('analysis.circulation.submonthly_analysis', []),
    "tracks": ("analysis.tracks.agency_comparison", []),
    "landfall": ("analysis.tracks.agency_comparison", ["landfall"]),
    "genesis": ("analysis.tracks.genesis_conditioned_decomposition", []),
    "speed": ("analysis.tracks.agency_comparison", ["speed"]),
    "steering": ("analysis.circulation.steering_layer_diagnostics", []),
    "trajectories": ("analysis.models.bam_primary_fit", []),
    "ocean": ("analysis.ocean.ohc_pipeline", []),
    "d26": ("analysis.ocean.d26_exposure", []),
    "pathways": ("analysis.tracks.pathway_closure", []),
    "figures": ("analysis.figures.manuscript_figures", []),
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
