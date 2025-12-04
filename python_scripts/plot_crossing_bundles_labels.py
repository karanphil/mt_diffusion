#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot a full (nb_bundles*nb_sections) x (nb_bundles*nb_sections)
group-averaged Dice matrix from multiple JSON files.
"""

import argparse
import logging
import json
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt

from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )

    p.add_argument('out_png',
                   help='Output PNG file for the Dice matrix plot.')

    p.add_argument('--in_jsons', nargs='+', required=True,
                   help='Input JSON files (one per subject).')

    p.add_argument('--nb_sections', type=int, default=20,
                   help='Number of sections in the bundle labels. '
                        'Default is 20.')
    
    p.add_argument('--bundles_names', nargs='+', default=None,
                   help='Subset of bundle names to include. '
                        'If not provided, all bundles are used.')
    
    p.add_argument('--highlight_threshold', type=float, default=0.2,
                   help='Threshold to highlight Dice scores in the matrix. '
                        'Default is 0.2.')

    p.add_argument('--save_txt', default=None,
                   help='Save the final averaged matrix as a .txt file. '
                        'This should be the path to the output .txt file.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


# Compute Dice matrix for ONE subject
def compute_subject_dice(crossing_info, bundle_names, nb_sections):

    nb_bundles = len(bundle_names)
    full_size = nb_bundles * nb_sections
    M = np.zeros((full_size, full_size))

    # ---- Preload voxel counts ----
    voxels = {}
    for b in bundle_names:
        voxels[b] = {}
        for s in range(1, nb_sections + 1):
            entries = crossing_info[b][str(s)]
            nb_vox = next((item[1] for item in entries if item[0] == "Nb_voxels"), None)
            if nb_vox is None:
                raise ValueError(f'"Nb_voxels" not found for {b}, section {s}')
            voxels[b][s] = nb_vox

    # ---- Fill Dice matrix ----
    for bi, bundle_i in enumerate(bundle_names):
        for bj, bundle_j in enumerate(bundle_names):

            row_offset = bi * nb_sections
            col_offset = bj * nb_sections

            for si in range(1, nb_sections + 1):
                entries = crossing_info[bundle_i][str(si)]

                for entry in entries:
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

    return M


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_jsons])
    assert_outputs_exist(parser, args, [args.out_png])

    # TODO : ne pas utiliser le DICE score!!!

    nb_sections = args.nb_sections
    dice_matrices = []
    bundle_names_ref = None

    # Load all subjects & compute Dice
    for json_file in args.in_jsons:
        with open(json_file, 'r') as f:
            crossing_info = json.load(f)

        all_bundles = list(crossing_info.keys())

        if args.bundles_names is not None:
            missing = set(args.bundles_names) - set(all_bundles)
            if missing:
                raise ValueError(f"Requested bundles not found in {json_file}: {missing}")
            bundle_names = args.bundles_names
        else:
            bundle_names = all_bundles

        # Enforce same ordering across subjects
        if bundle_names_ref is None:
            bundle_names_ref = bundle_names
        else:
            if bundle_names != bundle_names_ref:
                raise ValueError("Bundle order mismatch across subjects!")

        M = compute_subject_dice(crossing_info, bundle_names, nb_sections)
        dice_matrices.append(M)

    # Group mean Dice
    M = np.mean(dice_matrices, axis=0)

    # ---- Save matrix to TXT ----
    if args.save_txt is not None:
        np.savetxt(args.save_txt, M, fmt="%.6f")

    nb_bundles = len(bundle_names_ref)

    # Plot (UNCHANGED from your version)
    fig, ax = plt.subplots(figsize=(12, 12))
 
    im = ax.imshow(M, origin='lower', vmin=0, vmax=1,
                   extent=[-0.5, M.shape[1] - 0.5, -0.5, M.shape[0] - 0.5])

    # ---- Minor grid ----
    for k in range(1, nb_bundles * nb_sections):
        pos = k - 0.5
        ax.axhline(pos, linewidth=0.15, color='0.85', zorder=3)
        ax.axvline(pos, linewidth=0.15, color='0.85', zorder=3)

    # ---- Major grid ----
    for k in range(1, nb_bundles):
        pos = k * nb_sections - 0.5
        ax.axhline(pos, linewidth=0.6, color='0.3', zorder=3)
        ax.axvline(pos, linewidth=0.6, color='0.3', zorder=3)

    # ---- Highlight ----
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if M[i, j] >= args.highlight_threshold:
                rect = patches.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                         edgecolor='red', linewidth=0.2,
                                         zorder=10)
                ax.add_patch(rect)

    # ---- Bundle labels ----
    centers = np.arange(nb_bundles) * nb_sections + nb_sections / 2 - 0.5
    ax.set_xticks(centers)
    ax.set_yticks(centers)
    ax.set_xticklabels(bundle_names_ref, rotation=90)
    ax.set_yticklabels(bundle_names_ref)

    ax.set_xlabel("Bundles x Sections")
    ax.set_ylabel("Bundles x Sections")
    ax.set_title(f"Group mean Dice matrix (N={len(dice_matrices)})")

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, shrink=0.95)
    cbar.set_label("Dice score")

    plt.tight_layout()
    plt.savefig(args.out_png, dpi=1000)
    # plt.show()


if __name__ == "__main__":
    main()
