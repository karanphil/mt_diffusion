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

    p.add_argument('in_peak_values')

    p.add_argument('in_fixel_density_mask',
                   help='For one bundle.')

    p.add_argument('out_bundle_mtr')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_peak_values, args.in_fixel_density_mask])
    assert_outputs_exist(parser, args, [args.out_bundle_mtr])

    peak_values_img = nib.load(args.in_peak_values)
    peak_values = peak_values_img.get_fdata().astype(np.float32)

    fixel_density_img = nib.load(args.in_fixel_density_mask)
    fixel_density = fixel_density_img.get_fdata().astype(np.uint8)

    # Get the index of the first 1 along the last axis
    first_1_index = np.argmax(fixel_density, axis=3)

    # Create a mask where no 1s are found (i.e., all values were 0)
    no_1s_mask = ~np.any(fixel_density == 1, axis=3)

    # Prepare indices for advanced indexing
    X, Y, Z, _ = peak_values.shape
    i, j, k = np.indices((X, Y, Z))
    mtr = peak_values[i, j, k, first_1_index]

    # Optionally set mtr to 0 where no 1s were found
    mtr[no_1s_mask] = 0

    nib.save(nib.Nifti1Image(mtr, peak_values_img.affine), args.out_bundle_mtr)


if __name__ == "__main__":
    main()
