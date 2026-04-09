#!/usr/bin/env python3
"""
Compute the total radiation energy from mcscat and plot
the radial profile averaged over angles (theta, phi).

Usage:
    python mcscat_energy_vs_r.py [directory]
"""

import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import pandas as pd

# Parse through the directory
directory = sys.argv[1] if len(sys.argv) > 1 else '.'

if not os.path.isdir(directory):
    print(f"Error: {directory} is not a directory.")
    sys.exit(1)

# Find .athdf files
files = sorted(glob(os.path.join(directory, "*.athdf")))
if not files:
    print(f"No .athdf files found in {directory}")
    sys.exit(0)

print(f"Found {len(files)} .athdf files:")
for f in files:
    print("  ", os.path.basename(f))

# Process each file
for filename in files:
    print(f"\nProcessing {filename} ...")
    try:
        with h5py.File(filename, "r") as f:
            if "mcscat" not in f:
                print("  No 'mcscat' in file, skipping.")
                continue

            mc = np.array(f["mcscat"], dtype=np.float64)
            print(f"  mcscat shape = {mc.shape}")
            # (freq, r, theta, phi, angle)

            # Radial coordinate
            r = np.squeeze(np.array(f["x1v"], dtype=np.float64))
            #r = np.array(f["x1v"], dtype=np.float64)
            #print(np.array(f["x2v"], dtype=np.float64).shape)
            if r.shape[0] != mc.shape[1]:
                raise ValueError("x1v does not match radial dimension of mcscat")

        # ----------------------------------------------------------
        # Sum over frequency bins (axis 0) and angles (axis 4)
        # Average over theta (axis 2) and phi (axis 3)
        # ----------------------------------------------------------
        J_sum = np.sum(mc, axis=(4))          # -> (freq, r, theta, phi)
        
        J_sum = J_sum[0] #choose freq=0
        J_mean_ang = J_sum[:,0,0]
        #J_mean_ang = np.mean(J_sum, axis=(1,2)) # -> (r,)

        # Total integrated energy
        E_total = np.sum(mc)
        print(f"  Total mcscat energy (sum over all axes) = {E_total:.6e}")
        print(f" Mean mcscat value = {np.mean(mc):.3e}")
        # Plot <J>(r)
        plt.figure(figsize=(6,4))
        #if np.max(J_mean_ang) > 0:
        #    J_plot = J_mean_ang / np.max(J_mean_ang)
        #else:
        #    print("  Warning: mcscat is identically zero; skipping normalization.")
        J_plot = J_mean_ang
        r=r[:,0]
        #print(r[:,0])
        plt.plot(r, J_plot, marker='o',c='red',label='mcscat')
        plt.xlabel("r")
        plt.yscale('log')
        plt.ylabel(r"$J_{\theta=0,\phi=0}$")
        plt.title(f"J(r) from mcscat\n{os.path.basename(filename)}")
        plt.grid(True)

        df = pd.read_csv('./debugOutput/scatteringVals.csv',header=None,names=['x1','x2','x3','nu_mc','nu','scat'])
        
        df['r'] = np.sqrt(df['x1']**2 + df['x2']**2 + df['x3']**2)
        df = df[df['nu_mc'] == df['nu_mc'].min()] #choose one frequency bin
        plt.scatter(df['r'][::10]*1e11,df['scat'][::10],label="blacklight",marker='x')
        print(f"Mean scattering val as read by blacklight: {np.mean(df['scat']):.3e}")
        print(f"Mean scattering val near r=0: {np.mean(df['scat'][df['r']==df['r'].abs().min()]):.3e}")
        plt.legend()
        plt.tight_layout()

        #outname = os.path.splitext(filename)[0] + "_tegan_mcscat_vs_r.png"
        outname = "/PellaShared/kcu8rf/blacklight/plots/spherical_thomson/mscat_vs_r.png"
        plt.savefig(outname, dpi=150)
        plt.close()
        print(f"  Saved plot to {outname}")

    except Exception as e:
        print(f"  Error processing {filename}: {e}")
        continue

print("\nAll files processed successfully.")

