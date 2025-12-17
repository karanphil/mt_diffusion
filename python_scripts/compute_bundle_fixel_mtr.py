#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to compute the bundle-wise fixel-MTR for a single bundle.
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

    p.add_argument('in_fixel_mtr',
                   help='Fixel-MTR file, obtained from compute_fixel_mtr.py.')

    p.add_argument('in_fixel_density_mask',
                   help='Fixel density mask normalized by voxel for the '
                        'current bundle. This should be the result of '
                        'scil_bundle_fixel_analysis, named '
                        'as fixel_density_mask_voxel-norm_BUNDLE_NAME.nii.gz')
    
    p.add_argument('in_fixel_density_map',
                   help='Fixel density map normalized by voxel for the '
                        'current bundle. This should be the result of '
                        'scil_bundle_fixel_analysis, named '
                        'as fixel_density_map_voxel-norm_BUNDLE_NAME.nii.gz')

    p.add_argument('out_bundle_mtr',
                   help='Output bundle-wise fixel-MTR file.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_fixel_mtr,
                                 args.in_fixel_density_mask,
                                 args.in_fixel_density_map])
    assert_outputs_exist(parser, args, [args.out_bundle_mtr])

    peak_values_img = nib.load(args.in_fixel_mtr)
    peak_values = peak_values_img.get_fdata().astype(np.float32)

    fd_mask_img = nib.load(args.in_fixel_density_mask)
    fd_mask = fd_mask_img.get_fdata().astype(np.uint8)

    fd_map_img = nib.load(args.in_fixel_density_map)
    fd_map = fd_map_img.get_fdata().astype(np.float32)

    weights = fd_map * (fd_mask > 0)

    weight_sum = np.sum(weights, axis=3, keepdims=True)

    # Safe per-voxel normalization
    weights_norm = np.zeros_like(weights, dtype=np.float32)
    nonzero = weight_sum > 0

    print(peak_values.shape, weights.shape, weight_sum.shape, weights_norm.shape)

    weights_norm[nonzero] = weights[nonzero] / weight_sum[nonzero]

    mtr = np.sum(peak_values * weights_norm, axis=3)

    nib.save(nib.Nifti1Image(mtr, peak_values_img.affine), args.out_bundle_mtr)


if __name__ == "__main__":
    main()
