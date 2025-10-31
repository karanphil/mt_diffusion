#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to plot the track-profile of a bundle for both the MTR and the
fixel-wise MTR.
"""

import argparse
import logging

from cmcrameri import cm
from matplotlib.lines import Line2D 
import matplotlib as mpl
from matplotlib import pyplot as plt
import nibabel as nib
import numpy as np

from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('in_bundle_mtr')

    p.add_argument('in_bundle_fixel_mtr')

    p.add_argument('in_bundle_labels')

    p.add_argument('in_afd_fixel')

    p.add_argument('in_nufo')

    p.add_argument('out_dir')

    p.add_argument('--in_bundle_map')

    p.add_argument('--bundle_name')

    p.add_argument('--map_threshold', type=float, default=0.75,
                   help='Threshold to apply to the bundle map to create a '
                        'mask. Default is 0.1.')

    p.add_argument('--afd_threshold', type=float, default=0.2,
                   help='Threshold to apply to the AFD fixel map to discard '
                        'voxels with low AFD. Default is 0.2.')

    p.add_argument('--min_nvox', type=int, default=100,
                   help='Minimum number of voxels per bundle section to '
                        'consider it valid. Default is 100.')

    p.add_argument('--median', action='store_true',
                   help='Use the median instead of the mean to compute the '
                        'track-profile. Default is False.')
    
    p.add_argument('--variance', action='store_true',
                   help='Compute and display the variance of the MTR and '
                        'fixel-wise MTR for each bundle section. Default is '
                        'False.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_bundle_mtr, args.in_bundle_fixel_mtr,
                                 args.in_bundle_labels, args.in_afd_fixel,
                                 args.in_nufo])
    assert_outputs_exist(parser, args, [args.out_dir])

    mtr_img = nib.load(args.in_bundle_mtr)
    mtr = mtr_img.get_fdata().astype(np.float32)

    fixel_mtr_img = nib.load(args.in_bundle_fixel_mtr)
    fixel_mtr = fixel_mtr_img.get_fdata().astype(np.float32)

    labels_img = nib.load(args.in_bundle_labels)
    labels = labels_img.get_fdata().astype(np.uint8)

    afd_fixel_img = nib.load(args.in_afd_fixel)
    afd_fixel = afd_fixel_img.get_fdata().astype(np.float32)

    nufo_img = nib.load(args.in_nufo)
    nufo = nufo_img.get_fdata().astype(np.float32)

    if args.in_bundle_map:
        map_img = nib.load(args.in_bundle_map)
        map = map_img.get_fdata().astype(np.float32)
        mask = map >= args.map_threshold
    else:
        mask = np.ones(mtr.shape, dtype=bool)

    unique_labels = np.unique(labels[mask])
    unique_labels = unique_labels[unique_labels != 0]
    # unique_labels = unique_labels[unique_labels > 3]
    # unique_labels = unique_labels[unique_labels < 18]

    # cmap = ['orange', 'blue', 'red', 'green']
    cmap = cm.naviaS
    cmap_idx = np.arange(3, 1000, 1)
    norm = mpl.colors.Normalize(vmin=0.3, vmax=0.7)

    fixel_mtr_profile = np.zeros((len(unique_labels),))
    mtr_profile = np.zeros((len(unique_labels),))
    nufo_profile = np.zeros((len(unique_labels),))
    afd_profile = np.zeros((len(unique_labels),))
    for i, label in enumerate(unique_labels):
        label_mask = (labels == label) & mask & (afd_fixel > args.afd_threshold)
        print(label, np.sum(label_mask))
        if args.median:
            fixel_mtr_profile[i] = np.median(fixel_mtr[label_mask])
            mtr_profile[i] = np.median(mtr[label_mask])
            nufo_profile[i] = np.median(nufo[label_mask])
            afd_profile[i] = np.median(afd_fixel[label_mask])
        elif np.sum(afd_fixel[label_mask]) != 0 and np.sum(label_mask) >= args.min_nvox:
            fixel_mtr_profile[i] = np.average(fixel_mtr[label_mask],
                                              weights=afd_fixel[label_mask])
            mtr_profile[i] = np.average(mtr[label_mask],
                                        weights=afd_fixel[label_mask])
            nufo_profile[i] = np.average(nufo[label_mask],
                                         weights=afd_fixel[label_mask])
            afd_profile[i] = np.average(afd_fixel[label_mask])

    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(unique_labels[mtr_profile != 0], mtr_profile[mtr_profile != 0],
             label='MTR', marker='o', color=cmap(cmap_idx[0]))
    ax1.plot(unique_labels[fixel_mtr_profile != 0],
             fixel_mtr_profile[fixel_mtr_profile != 0], label='Fixel-wise MTR',
             marker='o', color=cmap(cmap_idx[1]))

    if args.variance:
        mtr_var = np.zeros((len(unique_labels),))
        fixel_mtr_var = np.zeros((len(unique_labels),))
        for i, label in enumerate(unique_labels):
            label_mask = (labels == label) & mask & (afd_fixel > args.afd_threshold)
            if np.sum(afd_fixel[label_mask]) != 0 and np.sum(label_mask) >= args.min_nvox:
                if args.median:
                    fixel_mtr_profile[i] = np.average(fixel_mtr[label_mask],
                                                    weights=afd_fixel[label_mask])
                    mtr_profile[i] = np.average(mtr[label_mask],
                                                weights=afd_fixel[label_mask])
                mtr_var[i] = np.average((mtr[label_mask]-mtr_profile[i])**2, 
                                        weights=afd_fixel[label_mask])
                fixel_mtr_var[i] = np.average((fixel_mtr[label_mask]-fixel_mtr_profile[i])**2,
                                            weights=afd_fixel[label_mask])
        ax1.fill_between(unique_labels[mtr_profile != 0],
                         mtr_profile[mtr_profile != 0] - np.sqrt(mtr_var[mtr_profile != 0]),
                         mtr_profile[mtr_profile != 0] + np.sqrt(mtr_var[mtr_profile != 0]),
                         color=cmap(cmap_idx[0]), alpha=0.2)
        ax1.fill_between(unique_labels[fixel_mtr_profile != 0],
                         fixel_mtr_profile[fixel_mtr_profile != 0] - np.sqrt(fixel_mtr_var[fixel_mtr_profile != 0]),
                         fixel_mtr_profile[fixel_mtr_profile != 0] + np.sqrt(fixel_mtr_var[fixel_mtr_profile != 0]),
                         color=cmap(cmap_idx[1]), alpha=0.2)

    # Add secondary axis for NuFO
    ax2 = ax1.twinx()
    ax2.plot(unique_labels[nufo_profile != 0],
             nufo_profile[nufo_profile != 0],
             linestyle='--', color="darkgrey", zorder=3)
    ax2.scatter(unique_labels[nufo_profile != 0],
                nufo_profile[nufo_profile != 0],
                c=afd_profile[nufo_profile != 0],
                cmap='Greys', norm=norm, edgecolors='darkgrey', zorder=4)
    complexity_handle = Line2D([0], [0],
                           color='darkgrey', linestyle='--',
                           marker='o', markerfacecolor='white',
                           markeredgecolor='darkgrey',
                           label='Complexity')
    ax2.set_ylabel('Mean NuFO', color="darkgrey")
    ax2.tick_params(axis='y', labelcolor="darkgrey")

    # Add AFD colorbar
    sm = plt.cm.ScalarMappable(cmap='Greys',
                               norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax2, pad=0.08)
    cbar.set_label('Mean AFD', color="darkgrey")
    cbar.ax.tick_params(color='darkgrey', labelcolor='darkgrey')
    cbar.outline.set_edgecolor('darkgrey')
    # cbar.tick_params(axis='y', labelcolor="darkgrey")

    # Axis labels and legend
    ax1.set_xlabel('Bundle section')
    ax1.set_ylabel('Mean MTR' if not args.median else 'Median MTR')
    bundle_name = f' for the {args.bundle_name} bundle' if args.bundle_name else ''
    ax1.set_title(f'Track-profile of MTR and fixel-wise MTR {bundle_name}')

    # Combine legends from both axes
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2 + [complexity_handle], labels_1 + labels_2 + ["Complexity"], loc='best')

    ax1.set_ylim(0.33, 0.45)
    ax1.set_xlim(0, 21)
    ax1.set_xticks(np.arange(1, 21, 1))
    ax2.set_ylim(1, 5)
    ax2.set_yticks(np.arange(1, 6, 1))

    plt.tight_layout()
    # plt.show()
    plt.savefig(args.out_dir + '/track_profile_{}.png'.format(args.bundle_name),
                dpi=300)


if __name__ == "__main__":
    main()
