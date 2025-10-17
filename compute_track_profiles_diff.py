#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to plot the absolute difference between test-retest of the track-profile
of a bundle for both the MTR and the fixel-wise MTR.
"""

import argparse
import logging

from cmcrameri import cm
from matplotlib import pyplot as plt
import nibabel as nib
import numpy as np

from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('in_bundle_mtr', nargs=2,
                   help='Input MTR bundle image (test and retest).')

    p.add_argument('in_bundle_fixel_mtr', nargs=2,
                   help='Input fixel-wise MTR bundle image (test and retest).')

    p.add_argument('in_bundle_labels', nargs=2)

    p.add_argument('in_afd_fixel', nargs=2)

    p.add_argument('out_dir')

    p.add_argument('--in_bundle_map', nargs=2)

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

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_bundle_mtr, args.in_bundle_fixel_mtr,
                                 args.in_bundle_labels, args.in_afd_fixel])
    assert_outputs_exist(parser, args, [args.out_dir])

    fixel_mtr_profile = np.zeros((2, 20))
    mtr_profile = np.zeros((2, 20))
    for j in range(2):
        mtr_img = nib.load(args.in_bundle_mtr[j])
        mtr = mtr_img.get_fdata().astype(np.float32)

        fixel_mtr_img = nib.load(args.in_bundle_fixel_mtr[j])
        fixel_mtr = fixel_mtr_img.get_fdata().astype(np.float32)

        labels_img = nib.load(args.in_bundle_labels[j])
        labels = labels_img.get_fdata().astype(np.uint8)

        afd_fixel_img = nib.load(args.in_afd_fixel[j])
        afd_fixel = afd_fixel_img.get_fdata().astype(np.float32)

        if args.in_bundle_map:
            map_img = nib.load(args.in_bundle_map[j])
            map = map_img.get_fdata().astype(np.float32)
            mask = map >= args.map_threshold
        else:
            mask = np.ones(mtr.shape, dtype=bool)

        unique_labels = np.unique(labels[mask])
        unique_labels = unique_labels[unique_labels != 0]

        cmap = cm.naviaS
        cmap_idx = np.arange(3, 1000, 1)

        for label in unique_labels:
            label_mask = (labels == label) & mask & (afd_fixel > args.afd_threshold)
            # print(label, np.sum(label_mask))
            if args.median:
                fixel_mtr_profile[j, label - 1] = np.median(fixel_mtr[label_mask])
                mtr_profile[j, label - 1] = np.median(mtr[label_mask])
            elif np.sum(afd_fixel[label_mask]) != 0 and np.sum(label_mask) >= args.min_nvox:
                fixel_mtr_profile[j, label - 1] = np.average(fixel_mtr[label_mask],
                                                weights=afd_fixel[label_mask])
                mtr_profile[j, label - 1] = np.average(mtr[label_mask],
                                            weights=afd_fixel[label_mask])
                
    diff_mtr_profile = np.abs(mtr_profile[0] - mtr_profile[1]) / mtr_profile[0] * 100
    diff_fixel_mtr_profile = np.abs(fixel_mtr_profile[0] - fixel_mtr_profile[1]) / fixel_mtr_profile[0] * 100

    fig, ax1 = plt.subplots(figsize=(8, 2))

    ax1.scatter(unique_labels[diff_mtr_profile > 0], diff_mtr_profile[diff_mtr_profile > 0],
             label='MTR', marker='o', color=cmap(cmap_idx[0]))
    ax1.scatter(unique_labels[diff_fixel_mtr_profile > 0],
             diff_fixel_mtr_profile[diff_fixel_mtr_profile > 0], label='Fixel-wise MTR',
             marker='o', color=cmap(cmap_idx[1]))

    # Axis labels and legend
    ax1.set_xlabel('Bundle section')
    ax1.set_ylabel('|%diff| scan-rescan')
    bundle_name = f' for the {args.bundle_name} bundle' if args.bundle_name else ''
    ax1.set_title(f'Track-profile repeatability {bundle_name}')

    # Combine legends from both axes
    # lines_1, labels_1 = ax1.get_legend_handles_labels()
    # ax1.legend(lines_1, labels_1, loc='best')

    ax1.set_ylim(0, 5)
    ax1.set_yticks(np.arange(1, 6, 1))
    ax1.set_xlim(0, 21)
    ax1.set_xticks(np.arange(1, 21, 1))

    plt.tight_layout()
    # plt.show()
    plt.savefig(args.out_dir + '/track_profile_{}_diff.png'.format(args.bundle_name),
                dpi=300)


if __name__ == "__main__":
    main()
