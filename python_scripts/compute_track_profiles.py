#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to compute the track-profile of a bundle for the MTR, fixel-wise MTR,
AFD and NuFO.
"""

import argparse
import logging

import nibabel as nib
import numpy as np

from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('in_bundle_mtr',
                   help='Input bundle-specific MTR image.')

    p.add_argument('in_bundle_fixel_mtr',
                   help='Input bundle-specific fixel-wise MTR image.')

    p.add_argument('in_bundle_labels',
                   help='Input bundle labels (sections) image.')

    p.add_argument('in_bundle_name', type=str,
                   help='Name of the bundle being processed.')

    p.add_argument('in_afd_fixel',
                   help='Input bundle-specific AFD fixel image.')

    p.add_argument('in_nufo',
                   help='Input NuFO image.')

    p.add_argument('out_dir',
                   help='Output directory to save the track-profiles in txt.')

    p.add_argument('--in_bundle_map',
                   help='Input bundle map to create a mask of the bundle. '
                        'If not provided, the whole bundle labels are used.')

    p.add_argument('--in_crossing_bundle_map',
                   help='Input crossing bundle map to create a mask of the '
                        'crossing bundle. This will produce track-profiles '
                        'using only the voxels comprised in this mask and the '
                        'usual masks. If not provided, the whole bundle '
                        'labels are used.')
    
    p.add_argument('--in_crossing_bundle_name',
                   help='Name of the crossing bundle.')

    p.add_argument('--map_threshold', type=float, default=1.0,
                   help='Threshold (higher-or-equal) to apply to the bundle '
                        'map to create a mask. Default is 1.0.')

    p.add_argument('--afd_threshold', type=float, default=0.3,
                   help='Threshold (higher-or-equal) to apply to the AFD '
                        'fixel map to discard voxels with low AFD. Default is '
                        '0.3.')

    p.add_argument('--min_nvox', type=int, default=100,
                   help='Minimum number of voxels per bundle section to '
                        'consider it valid. Default is 100.')
    
    p.add_argument('--nb_sections', type=int, default=20,
                   help='Number of sections in the bundle labels. '
                        'Default is 20.')

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
    afd_mask = afd_fixel >= args.afd_threshold

    nufo_img = nib.load(args.in_nufo)
    nufo = nufo_img.get_fdata().astype(np.float32)

    if args.in_bundle_map:
        map_img = nib.load(args.in_bundle_map)
        map = map_img.get_fdata().astype(np.float32)
        mask = map >= args.map_threshold
    else:
        mask = np.ones(mtr.shape, dtype=bool)

    if args.in_crossing_bundle_map:
        map_img = nib.load(args.in_crossing_bundle_map)
        map = map_img.get_fdata().astype(np.float32)
        crossing_mask = map >= args.map_threshold
    else:
        crossing_mask = np.ones(mtr.shape, dtype=bool)

    labels_list = np.arange(1, args.nb_sections + 1, 1)

    fixel_mtr_profile = np.zeros((len(labels_list),))
    mtr_profile = np.zeros((len(labels_list),))
    nufo_profile = np.zeros((len(labels_list),))
    afd_profile = np.zeros((len(labels_list),))
    nb_voxels_profile = np.zeros((len(labels_list),))
    for i, label in enumerate(labels_list):
        label_mask = (labels == label) & mask & afd_mask & crossing_mask
        nb_voxels_profile[i] = np.sum(label_mask)
        if np.sum(afd_fixel[label_mask]) != 0 and nb_voxels_profile[i] >= args.min_nvox:
            fixel_mtr_profile[i] = np.average(fixel_mtr[label_mask],
                                              weights=afd_fixel[label_mask])
            mtr_profile[i] = np.average(mtr[label_mask],
                                        weights=afd_fixel[label_mask])
            nufo_profile[i] = np.average(nufo[label_mask],
                                         weights=afd_fixel[label_mask])
            afd_profile[i] = np.average(afd_fixel[label_mask])

    if args.in_crossing_bundle_map and args.in_crossing_bundle_name:
        ext_name = f"_crossing_{args.in_crossing_bundle_name}"
    elif args.in_crossing_bundle_map:
        ext_name = "_crossing"
    else:
        ext_name = ""

    np.savetxt(f"{args.out_dir}/mtr_profile_{args.in_bundle_name}{ext_name}.txt",
               mtr_profile)
    np.savetxt(f"{args.out_dir}/fixel_mtr_profile_{args.in_bundle_name}{ext_name}.txt",
               fixel_mtr_profile)
    np.savetxt(f"{args.out_dir}/nufo_profile_{args.in_bundle_name}{ext_name}.txt",
               nufo_profile)
    np.savetxt(f"{args.out_dir}/afd_profile_{args.in_bundle_name}{ext_name}.txt",
               afd_profile)
    np.savetxt(f"{args.out_dir}/nb_voxels_profile_{args.in_bundle_name}{ext_name}.txt",
               nb_voxels_profile)


if __name__ == "__main__":
    main()
