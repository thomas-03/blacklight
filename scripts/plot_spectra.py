#!/usr/bin/env python3
"""
plot_spectra.py

Command-line-driven version of the original plot_spec.py: plots a spectrum
from a single blacklight .npz image, optionally overlaid with a Monte Carlo
spectrum (bap.MCSpec) and/or a blackbody comparison curve.

Intended to be called once per test case by run_tests.sh, e.g.:

  Non-MC case:
    python3 plot_spectra.py \
        --image ./output/example_isoth.npz \
        --output ./output/example_isoth_spec.png \
        --label example_isoth

  MC case (mc-dir comes from the "mc_file" line in the .input file):
    python3 plot_spectra.py \
        --image ./output/example_mc_isoth_gr.npz \
        --output ./output/example_mc_isoth_gr_spec.png \
        --label example_mc_isoth_gr \
        --mc-dir /PellaShared/kcu8rf/sks_ring_comp_lowphot/
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless-safe for batch/script runs
import matplotlib.pyplot as plt
import bap


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot a spectrum from a blacklight .npz image, "
                     "optionally overlaid with an MC spectrum and/or blackbody curve."
    )
    p.add_argument("--image", required=True,
                   help="Path to the blacklight .npz image file (bap.Image input).")
    p.add_argument("--output", required=True,
                   help="Path to save the output plot, e.g. spectrum.png")
    p.add_argument("--label", default=None,
                   help="Legend label for the blacklight curve. Defaults to the image filename.")

    # --- Monte Carlo spectrum options (only used when --mc-dir is passed) ---
    p.add_argument("--mc-dir", default=None,
                   help="Directory of MC run outputs, passed to bap.MCSpec. "
                        "Only pass this for example_mc_* test cases.")
    p.add_argument("--mc-nproc", type=int, default=16,
                   help="Number of processes for bap.MCSpec (default: 16).")
    p.add_argument("--mc-screen", default="inner_radius",
                   help="`screen` argument for bap.MCSpec (default: inner_radius).")
    p.add_argument("--mc-nfreq", type=int, default=40,
                   help="Number of frequency bins for bap.MCSpec (default: 40).")
    p.add_argument("--mc-emin", type=float, default=1,
                   help="Minimum energy for bap.MCSpec (default: 1).")
    p.add_argument("--mc-emax", type=float, default=10000,
                   help="Maximum energy for bap.MCSpec (default: 10000).")
    p.add_argument("--mc-overwrite", dest="mc_overwrite", action="store_true", default=False,
                   help="Overwrite a cached MC spectrum (default is to not overwrite).")

    # --- Blackbody comparison curve ---
    p.add_argument("--no-blackbody", dest="plot_blackbody", action="store_false", default=True,
                   help="Disable the blackbody comparison curve (enabled by default).")
    p.add_argument("--temperature", type=float, default=1e7,
                   help="Blackbody temperature in K (default: 1e7).")
    p.add_argument("--radius", type=float, default=8859750.2283,
                   help="Emitting radius in cm, used as area = 2*pi*radius**2 (default matches "
                        "the original plot_spec.py value).")
    p.add_argument("--redshift_re", type=float, default=0.0,
                   help="Gravitational Radii for redshift for the blackbody curve (default: 6.0).")

    # --- Axis options ---
    p.add_argument("--ylim", type=float, nargs=2, default=[1e30, 1e40],
                   help="Y-axis limits, e.g. --ylim 1e34 1e38 (default: 1e34 1e38).")
    p.add_argument("--xlim", type=float, nargs=2, default=None,
                   help="Optional x-axis limits, e.g. --xlim 1e1 1e4.")

    return p.parse_args()


def main():
    args = parse_args()

    im = bap.Image(args.image)
    fig, ax = plt.subplots()

    mc = None
    if args.mc_dir:
        mc = bap.MCSpec(
            args.mc_dir,
            args.mc_nproc,
            screen=args.mc_screen,
            nfreq=args.mc_nfreq,
            emin=args.mc_emin,
            emax=args.mc_emax,
            overwrite=args.mc_overwrite,
        )

    area = 2 * np.pi * (args.radius ** 2)
    label = args.label or args.image

    bap.plot_spectra(
        [im],
        ax,
        MC_spec=mc,
        plot_blackbody=args.plot_blackbody,
        temperature=args.temperature,
        area=area,
        labels=[label],
    )

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylim(*args.ylim)
    if args.xlim:
        ax.set_xlim(*args.xlim)
    ax.grid()

    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Saved plot: {args.output}")


if __name__ == "__main__":
    main()