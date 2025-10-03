#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to plot the track-profile of a bundle for both the MTR and the
fixel-wise MTR.
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

    p.add_argument('in_bundle_mtr')

    p.add_argument('in_bundle_fixel_mtr')

    p.add_argument('in_bundle_labels')

    p.add_argument('out_dir')

    p.add_argument('--in_bundle_map')

    p.add_argument('--bundle_name')

    p.add_argument('--map_threshold', type=float, default=0.75,
                   help='Threshold to apply to the bundle map to create a '
                        'mask. Default is 0.1.')
    
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
    unique_labels = unique_labels[unique_labels > 3]
    unique_labels = unique_labels[unique_labels < 18]

    # cmap = ['orange', 'blue', 'red', 'green']
    cmap = cm.naviaS
    cmap_idx = np.arange(3, 1000, 1)

    fixel_mtr_profile = np.zeros((len(unique_labels),))
    mtr_profile = np.zeros((len(unique_labels),))
    for i, label in enumerate(unique_labels):
        label_mask = (labels == label) & mask
        print(label, np.sum(label_mask))
        if args.median:
            fixel_mtr_profile[i] = np.median(fixel_mtr[label_mask])
            mtr_profile[i] = np.median(mtr[label_mask])
        else:
            fixel_mtr_profile[i] = np.mean(fixel_mtr[label_mask])
            mtr_profile[i] = np.mean(mtr[label_mask])

    plt.plot(unique_labels, mtr_profile, label='MTR', marker='o',
            color=cmap(cmap_idx[0]))
    plt.plot(unique_labels, fixel_mtr_profile, label='Fixel-wise MTR',
            marker='o', color=cmap(cmap_idx[1]))

    if args.variance:
        mtr_var = np.zeros((len(unique_labels),))
        fixel_mtr_var = np.zeros((len(unique_labels),))
        for i, label in enumerate(unique_labels):
            label_mask = (labels == label) & mask
            mtr_var[i] = np.var(mtr[label_mask])
            fixel_mtr_var[i] = np.var(fixel_mtr[label_mask])
            if args.median:
                fixel_mtr_profile[i] = np.mean(fixel_mtr[label_mask])
                mtr_profile[i] = np.mean(mtr[label_mask])
        plt.fill_between(unique_labels,
                         mtr_profile - np.sqrt(mtr_var),
                         mtr_profile + np.sqrt(mtr_var),
                         color=cmap(cmap_idx[0]), alpha=0.2)
        plt.fill_between(unique_labels,
                         fixel_mtr_profile - np.sqrt(fixel_mtr_var),
                         fixel_mtr_profile + np.sqrt(fixel_mtr_var),
                         color=cmap(cmap_idx[1]), alpha=0.2)

    plt.xlabel('Bundle section')
    if args.median:
        plt.ylabel('Median MTR')
    else:
        plt.ylabel('Mean MTR')
    bundle_name = ' for the {} bundle'.format(args.bundle_name)
    plt.title('Track-profile of MTR and fixel-wise MTR' + bundle_name)
    plt.legend()
    plt.ylim(0.33, 0.45)
    plt.xlim(3, 18)
    # plt.show()
    plt.savefig(args.out_dir + '/track_profile_{}.png'.format(args.bundle_name),
                dpi=300)


if __name__ == "__main__":
    main()
