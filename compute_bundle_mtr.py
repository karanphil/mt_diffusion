#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

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

    p.add_argument('in_mtr')

    p.add_argument('in_voxel_density_masks')

    p.add_argument('in_LUT')

    p.add_argument('out_dir')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_peak_values, args.in_fixel_density_mask])
    assert_outputs_exist(parser, args, [args.out_dir])

    mtr_img = nib.load(args.in_peak_values)
    mtr = mtr_img.get_fdata().astype(np.float32)

    voxel_density_img = nib.load(args.in_fixel_density_mask)
    bundles_mask = voxel_density_img.get_fdata().astype(bool)

    lut = np.loadtxt(args.in_LUT, dtype=str)

    for i in range(bundles_mask.shape[-1]):
        bundle_mtr = np.zeros(mtr.shape)
        bundle_mtr[bundles_mask[..., i]] = mtr[bundles_mask[..., i]]

        bundle_name = lut[0][i]

        nib.save(nib.Nifti1Image(bundle_mtr, mtr_img.affine), args.out_dir + f'/mtr_{bundle_name}.nii.gz')


if __name__ == "__main__":
    main()
