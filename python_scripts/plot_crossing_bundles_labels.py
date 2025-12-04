#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot a full (nb_bundles*nb_sections) x (nb_bundles*nb_sections)
matrix where each block (20x20) represents bundle crossings.
Shared voxels are normalized by min(Ni, Nj) by default.
"""

import argparse
import json
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )

    p.add_argument('in_json',
                   help='Input JSON file produced by the crossing script.')

    p.add_argument('--nb_sections', type=int, default=20,
                   help='Number of sections in the bundle labels. '
                        'Default is 20.')

    p.add_argument('--vmax', type=float, default=1.0,
                   help='Upper bound for colormap (default: 1.0).')
    
    # TODO: add option to select a subset of bundles to plot

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    # Load JSON
    with open(args.in_json, 'r') as f:
        crossing_info = json.load(f)

    bundle_names = list(crossing_info.keys())
    nb_bundles = len(bundle_names)
    nb_sections = args.nb_sections

    full_size = nb_bundles * nb_sections
    M = np.zeros((full_size, full_size))

    # Preload voxel counts USING "Nb_voxels" KEY
    voxels = {}

    for b in bundle_names:
        voxels[b] = {}
        for s in range(1, nb_sections + 1):
            entries = crossing_info[b][str(s)]

            # Explicitly find ("Nb_voxels", value)
            nb_vox = None
            for item in entries:
                if item[0] == "Nb_voxels":
                    nb_vox = item[1]
                    break

            if nb_vox is None:
                raise ValueError(f'"Nb_voxels" not found for {b}, section {s}')

            voxels[b][s] = nb_vox

    # Fill global matrix
    for bi, bundle_i in enumerate(bundle_names):
        for bj, bundle_j in enumerate(bundle_names):

            row_offset = bi * nb_sections
            col_offset = bj * nb_sections

            for si in range(1, nb_sections + 1):
                entries = crossing_info[bundle_i][str(si)]

                for entry in entries:
                    # Skip the "Nb_voxels" entry
                    if entry[0] == "Nb_voxels":
                        continue

                    bname, sj, shared = entry

                    if bname != bundle_j:
                        continue

                    Ni = voxels[bundle_i][si]
                    Nj = voxels[bundle_j][sj]

                    dice = 2.0 * shared / (Ni + Nj) if (Ni + Nj) > 0 else 0.0

                    ii = row_offset + (si - 1)
                    jj = col_offset + (sj - 1)

                    M[ii, jj] = dice

    # Plot
    fig, ax = plt.subplots(figsize=(12, 12))

    im = ax.imshow(M, origin='lower', vmin=0, vmax=args.vmax, extent=[-0.5, M.shape[1] - 0.5, -0.5, M.shape[0] - 0.5])

    # ---- Minor grid: section blocks ----
    for k in range(1, nb_bundles * nb_sections):
        pos = k - 0.5
        ax.axhline(pos, linewidth=0.15, color='0.85', zorder=3)
        ax.axvline(pos, linewidth=0.15, color='0.85', zorder=3)

    # ---- Major grid: bundle blocks ----
    for k in range(1, nb_bundles):
        pos = k * nb_sections - 0.5
        ax.axhline(pos, linewidth=0.6, color='0.3', zorder=3)
        ax.axvline(pos, linewidth=0.6, color='0.3', zorder=3)

    # ---- Highlight cells with dice score above threshold ----
    threshold = 0.2
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if M[i, j] >= threshold:
                rect = patches.Rectangle(
                    (j - 0.5, i - 0.5),  # lower-left corner
                    1, 1,               # width, height = 1 cell
                    fill=False,
                    edgecolor='red',    # highlight color
                    linewidth=0.2,
                    zorder=10
                )
                ax.add_patch(rect)

    # ---- Bundle labels at block centers ----
    centers = np.arange(nb_bundles) * nb_sections + nb_sections / 2 - 0.5
    ax.set_xticks(centers)
    ax.set_yticks(centers)
    ax.set_xticklabels(bundle_names, rotation=90)
    ax.set_yticklabels(bundle_names)

    # ---- Axis labels ----
    ax.set_xlabel("Bundles x Sections")
    ax.set_ylabel("Bundles x Sections")
    ax.set_title("Global shared-voxel dice score matrix")

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, shrink=0.95)
    cbar.set_label("Dice score")

    plt.tight_layout()
    #plt.show()
    plt.savefig('crossing_bundles_matrix.png', dpi=1000)


if __name__ == "__main__":
    main()
