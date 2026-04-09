"""
compare_J_vs_Er.py

Written/readjusted script for the spherical run
current status: not good, being silly, needs to be adjusted.

Script Objective:
Compute the total radiation energy from mcscat and plot
the radial profile averaged over angles (theta, phi).
"""

import sys
import os
import h5py
import athena_read
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

directory = sys.argv[1] if len(sys.argv) > 1 else '.'

files = sorted(glob(os.path.join(directory, "*.athdf")))
if not files:
    print("No .athdf files found.")
    sys.exit(0)

for filename in files:
    print(f"\nProcessing {filename} ...")

    """
    with h5py.File(filename, "r") as f:
    if "mcscat" not in f:
        print("  No mcscat field, skipping.")
        continue
    """
    mcdict = athena_read.athdf(filename,raw=True)
    print(mcdict.keys())
    mcscat = mcdict["sscat0"]
    r = mcdict["x1v"]
    print(mcscat)
    J_plot = mcscat[0,0,:]
 
        
    """    
            mc = np.array(f["mcscat"], dtype=np.float64)
            print("  mcscat shape =", mc.shape)
            # Expected: (nfreq, nr, ntheta, nphi, nangle)

            # Robust radial coordinate read
            r_raw = np.array(f["x1v"], dtype=np.float64)
            r = np.squeeze(r_raw)

            # Sanity check, wahoowa
            if r.shape[0] != mc.shape[1]:
                raise ValueError(
                    f"x1v length {r.shape[0]} does not match mcscat radial size {mc.shape[1]}"
                )

        # Sum over both the frequency and the angle
        # review this
        J_sum = np.sum(mc, axis=(0, 4))         # -> (r, theta, phi)
        J_mean_ang = np.mean(J_sum, axis=(1,2))  # -> (r,)

        E_total = np.sum(mc)
        print(f"  Total mcscat energy = {E_total:.6e}")

        if np.max(J_mean_ang) > 0:
            J_plot = J_mean_ang / np.max(J_mean_ang)
        else:
            print("  Warning: mcscat is identically zero; skipping normalization.")
            J_plot = J_mean_ang
    """
    plt.figure(figsize=(6,4))
    plt.plot(r, J_plot, marker='o')
    plt.xlabel("r")
    plt.ylabel(r"Normalized $\langle J \rangle_{\theta,\phi}$")
    plt.title(os.path.basename(filename))
    plt.grid(True)
    plt.tight_layout()

    outname = os.path.splitext(filename)[0] + "_mcscat_vs_r.png"
    plt.savefig(outname, dpi=150)
    plt.close()
    print(f"  Saved plot to {outname}")
"""
except Exception as e:
    print(f"  Error processing {filename}: {e}")
"""

