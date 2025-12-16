#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot a full (nb_bundles*nb_sections) x (nb_bundles*nb_sections)
group-averaged overlap matrix from multiple JSON files.
"""

import argparse
import logging
import json
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt

from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


# ---- Bundle family exclusions ----
EXCLUDE_FAMILIES = [
    ["SLF"],      # all SLF parts exclude each other
    ["CC"],       # all CC parts exclude each other
]

# ---- Specific pair exclusions (parallel tracts) ----
EXCLUDE_PAIRS = [
    ("AF", "SLF"), ("CG", "SLF"), ("ILF", "OR")
]


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )

    p.add_argument('out_png',
                   help='Output PNG file for the overlap matrix plot.')

    p.add_argument('--in_jsons', nargs='+', required=True,
                   help='Input JSON files (one per subject).')

    p.add_argument('--nb_sections', type=int, default=20,
                   help='Number of sections in the bundle labels. '
                        'Default is 20.')
    
    p.add_argument('--bundles_names', nargs='+', default=None,
                   help='Subset of bundle names to include. '
                        'If not provided, all bundles are used.')
    
    # p.add_argument('--highlight_threshold', type=float, default=50,
    #                help='Threshold to highlight percentage overlap in the '
    #                     'matrix. Default is 50%.')

    p.add_argument('--min_nvox', type=int, default=100,
                   help='Minimum number of voxels per bundle section to '
                        'consider it valid. Default is 100.')

    p.add_argument('--save_full_txt', default=None,
                   help='Save the final averaged matrix as a .txt file. '
                        'This should be the path to the output .txt file.')
    
    p.add_argument('--save_crossing_txt', default=None,
                   help='Save all section-to-bundle crossings no matter the '
                        'percentage. This should be the path to the output '
                        '.txt file.')
    
    p.add_argument('--save_important_txt', default=None,
                   help='Save a .txt file listing section-to-bundle '
                        'crossings greater than 10%, excluding same-family '
                        'and parallel bundles. This should be the path to the '
                        'output .txt file.')
    
    p.add_argument('--important_threshold', type=float, default=20.0,
                   help='Threshold for important crossings to be saved '
                        'in the text file. Default is 20%.')

    p.add_argument('--use_crameri', action='store_true',
                   help='Use the Crameri colormap for plotting. '
                        'Default is False.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


# Compute overlap matrix for ONE subject
def compute_subject_overlap(crossing_info, bundle_names, nb_sections, min_nvox):
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

    # ---- Fill overlap matrix ----
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
                    # Nj = voxels[bundle_j][sj]
                    overlap = shared / Ni * 100 if Ni >= min_nvox else np.nan
                    # dice = 2.0 * shared / (Ni + Nj) if (Ni + Nj) > 0 else 0.0
                    ii = row_offset + (si - 1)
                    jj = col_offset + (sj - 1)
                    M[ii, jj] = overlap

    return M


def is_excluded(bundle_a, bundle_b):
    # Same family exclusion
    for fam in EXCLUDE_FAMILIES:
        if any(f in bundle_a for f in fam) and any(f in bundle_b for f in fam):
            return True

    # Explicit pair exclusion (both directions)
    for a, b in EXCLUDE_PAIRS:
        if (a in bundle_a and b in bundle_b) or (b in bundle_a and a in bundle_b):
            return True

    return False


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_jsons])
    assert_outputs_exist(parser, args, [args.out_png])

    nb_sections = args.nb_sections
    overlap_matrices = []
    bundle_names_ref = None

    # Load all subjects & compute overlap
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

        M = compute_subject_overlap(crossing_info, bundle_names, nb_sections,
                                    args.min_nvox)
        overlap_matrices.append(M)

    # Group mean overlap
    M = np.nanmean(overlap_matrices, axis=0)
    M = np.nan_to_num(M, nan=0.0)

    # Compute section-to-bundle crossing summary
    nb_bundles = len(bundle_names_ref)
    C = np.zeros((nb_bundles * nb_sections, nb_bundles))
    for bi in range(nb_bundles):
        for si in range(nb_sections):
            row = bi * nb_sections + si
            for bj in range(nb_bundles):
                col_start = bj * nb_sections
                col_end = col_start + nb_sections
                # Sum over the 20 sections of bundle bj
                C[row, bj] = np.sum(M[row, col_start:col_end])

    # Save matrix to TXT
    if args.save_full_txt is not None:
        np.savetxt(args.save_full_txt, M, fmt="%.6f")

    # Save crossings > args.important_threshold to TXT
    if args.save_important_txt is not None:
        C_threshold = args.important_threshold
        with open(args.save_important_txt, "w") as f:
            f.write("Section-to-bundle crossings > {}% (group-averaged)\n".format(args.important_threshold))
            f.write("Excluded: same-family bundles + parallel bundles\n\n")
            for bi, bundle_i in enumerate(bundle_names_ref):
                for si in range(nb_sections):
                    row = bi * nb_sections + si
                    source_label = f"{bundle_i} - section {si+1}"
                    for bj, bundle_j in enumerate(bundle_names_ref):
                        # ---- Exclusions ----
                        if is_excluded(bundle_i, bundle_j):
                            continue
                        val = C[row, bj]
                        if val >= C_threshold:
                            f.write(f"{source_label:<35} --> {bundle_j:<20} : {val:.2f}%\n")
    # Save all crossings to TXT
    if args.save_crossing_txt is not None:
        with open(args.save_crossing_txt, "w") as f:
            f.write("All section-to-bundle crossings (group-averaged)\n\n")
            f.write("Excluded: same-family bundles + parallel bundles\n\n")
            for bi, bundle_i in enumerate(bundle_names_ref):
                for si in range(nb_sections):
                    row = bi * nb_sections + si
                    source_label = f"{bundle_i} - section {si+1}"
                    for bj, bundle_j in enumerate(bundle_names_ref):
                        # ---- Exclusions ----
                        if is_excluded(bundle_i, bundle_j):
                            continue
                        val = C[row, bj]
                        if val >= 1:
                            f.write(f"{source_label:<35} --> {bundle_j:<20} : {val:.2f}%\n")

    # Choose colormap
    if args.use_crameri:
        from cmcrameri import cm
        cmap = cm.navia
    else:
        cmap = plt.get_cmap('bone')

    # Plot section-to-bundle matrix
    fig1, ax1 = plt.subplots(figsize=(12, 6))
    im1 = ax1.imshow(C.T, origin='lower', aspect='auto', cmap=cmap, vmin=0,
                     vmax=100)

    # ---- X axis = bundles × sections ----
    centers_x = np.arange(nb_bundles) * nb_sections + nb_sections / 2 - 0.5
    ax1.set_xticks(centers_x)
    ax1.set_xticklabels(bundle_names_ref, rotation=90)

    # ---- Y axis = bundles ----
    ax1.set_yticks(np.arange(nb_bundles))
    ax1.set_yticklabels(bundle_names_ref)

    # ---- Vertical major grid between bundles ----
    for k in range(1, nb_bundles):
        pos = k * nb_sections - 0.5
        ax1.axvline(pos, linewidth=0.6, color='0.3')

    # ---- Horizontal major grid between bundles ----
    for k in range(1, nb_bundles):
        ax1.axhline(k - 0.5, linewidth=0.6, color='0.3')

    ax1.set_xlabel("Source bundle (split by sections)")
    ax1.set_ylabel("Target bundle (averaged over sections)")
    # ax1.set_title("Bundle-by-bundle cumulative section overlap")

    cbar1 = fig1.colorbar(im1, ax=ax1, fraction=0.05, pad=0.02)
    cbar1.set_label("Percentage overlap")

    plt.tight_layout()
    plt.savefig(args.out_png, dpi=1000)
    # plt.show()

    # Plot section-to-section matrix
    fig2, ax2 = plt.subplots(figsize=(12, 12))
    im2 = ax2.imshow(M.T, cmap=cmap, origin='lower', vmin=0, vmax=100,
                     extent=[-0.5, M.shape[1] - 0.5, -0.5, M.shape[0] - 0.5])

    # ---- Minor grid ----
    for k in range(1, nb_bundles * nb_sections):
        pos = k - 0.5
        ax2.axhline(pos, linewidth=0.05, color='0.3', zorder=3)
        ax2.axvline(pos, linewidth=0.05, color='0.3', zorder=3)

    # ---- Major grid ----
    for k in range(1, nb_bundles):
        pos = k * nb_sections - 0.5
        ax2.axhline(pos, linewidth=0.6, color='0.3', zorder=3)
        ax2.axvline(pos, linewidth=0.6, color='0.3', zorder=3)

    # # ---- Highlight ----
    # for i in range(M.T.shape[0]):
    #     for j in range(M.T.shape[1]):
    #         if M.T[i, j] >= args.highlight_threshold:
    #             rect = patches.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
    #                                      edgecolor='red', linewidth=0.1,
    #                                      zorder=10)
    #             ax.add_patch(rect)

    # ---- Bundle labels ----
    centers = np.arange(nb_bundles) * nb_sections + nb_sections / 2 - 0.5
    ax2.set_xticks(centers)
    ax2.set_yticks(centers)
    ax2.set_xticklabels(bundle_names_ref, rotation=90)
    ax2.set_yticklabels(bundle_names_ref)

    ax2.set_xlabel("Source bundle (split by sections)")
    ax2.set_ylabel("Target bundle (split by sections)")
    # ax2.set_title(f"Group mean overlap matrix")

    cbar2 = fig2.colorbar(im2, ax=ax2, fraction=0.035, pad=0.02, shrink=0.95)
    cbar2.set_label("Percentage overlap")

    plt.tight_layout()
    plt.savefig(args.out_png.replace(".png", "_full.png"), dpi=1000)
    # plt.show()


if __name__ == "__main__":
    main()
