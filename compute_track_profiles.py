#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to plot the track-profile of a bundle for both the MTR and the
fixel-wise MTR.
"""

import argparse
import logging

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

    p.add_argument('out_dir')

    p.add_argument('--in_bundle_map')

    p.add_argument('--map_threshold', type=float, default=0.75,
                   help='Threshold to apply to the bundle map to create a '
                        'mask. Default is 0.1.')
    
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
                                 args.in_bundle_labels])
    assert_outputs_exist(parser, args, [args.out_dir])

    mtr_img = nib.load(args.in_bundle_mtr)
    mtr = mtr_img.get_fdata().astype(np.float32)

    fixel_mtr_img = nib.load(args.in_bundle_fixel_mtr)
    fixel_mtr = fixel_mtr_img.get_fdata().astype(np.float32)

    labels_img = nib.load(args.in_bundle_labels)
    labels = labels_img.get_fdata().astype(np.uint8)

    if args.in_bundle_map:
        map_img = nib.load(args.in_bundle_map)
        map = map_img.get_fdata().astype(np.float32)
        mask = map >= args.map_threshold
    else:
        mask = np.ones(mtr.shape, dtype=bool)

    unique_labels = np.unique(labels[mask])
    unique_labels = unique_labels[unique_labels != 0]

    fixel_mtr_profile = np.zeros((len(unique_labels),))
    mtr_profile = np.zeros((len(unique_labels),))
    for i, label in enumerate(unique_labels):
        label_mask = (labels == label) & mask
        if args.median:
            fixel_mtr_profile[i] = np.median(fixel_mtr[label_mask])
            mtr_profile[i] = np.median(mtr[label_mask])
        else:
            fixel_mtr_profile[i] = np.mean(fixel_mtr[label_mask])
            mtr_profile[i] = np.mean(mtr[label_mask])
        

    plt.plot(unique_labels, mtr_profile, label='MTR', marker='o',
             color='orange')
    plt.plot(unique_labels, fixel_mtr_profile, label='Fixel-wise MTR',
             marker='^', color='blue')
    plt.xlabel('Bundle section')
    if args.median:
        plt.ylabel('Median MTR')
    else:
        plt.ylabel('Mean MTR')
    plt.title('Track-profile of MTR and fixel-wise MTR')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
