#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to compute the bundle-wise MTR for all bundles.
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

    p.add_argument('in_mtr',
                   help='Input MTR file.')

    p.add_argument('in_voxel_density_masks',
                   help='Voxel density masks normalized by voxel for '
                        'each bundle. This should be the result of '
                        'scil_bundle_fixel_analysis, named '
                        'as voxel_density_masks_voxel-norm.nii.gz')

    p.add_argument('in_LUT',
                   help='Look-Up Table (LUT) file outputed by '
                        'scil_bundle_fixel_analysis, named '
                        'as bundles_LUT.txt')

    p.add_argument('out_dir',
                   help='Output directory to save the bundle-wise MTR files.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_mtr, args.in_voxel_density_masks])
    assert_outputs_exist(parser, args, [args.out_dir])

    mtr_img = nib.load(args.in_mtr)
    mtr = mtr_img.get_fdata().astype(np.float32)

    voxel_density_img = nib.load(args.in_voxel_density_masks)
    bundles_mask = voxel_density_img.get_fdata().astype(bool)

    lut = np.loadtxt(args.in_LUT, dtype=str)

    for i in range(bundles_mask.shape[-1]):
        bundle_mtr = np.zeros(mtr.shape)
        bundle_mtr[bundles_mask[..., i]] = mtr[bundles_mask[..., i]]

        bundle_name = lut[0][i]

        nib.save(nib.Nifti1Image(bundle_mtr, mtr_img.affine),
                 args.out_dir + f'/mtr_{bundle_name}.nii.gz')


if __name__ == "__main__":
    main()
