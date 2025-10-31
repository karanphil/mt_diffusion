#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to compare the fixel-MTR values of the two first peaks in a voxel.
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

    p.add_argument('out_differences',
                   help='Output image containing the difference between '
                        'the two first peaks MTR values.')

    p.add_argument('out_diff_mask',
                   help='Output binary mask where the difference between '
                        'the two first peaks is significant.')

    p.add_argument('out_no_diff_mask',
                    help='Output binary mask where there is no significant '
                         'difference between the two first peaks.')

    p.add_argument('out_crossing_mask',
                   help='Output binary mask where there are at least two '
                        'peaks with MTR values above the minimum MTR.')

    p.add_argument('--min_mtr', type=float, default=0.2,
                   help='Minimum MTR value to consider a peak valid.')
    
    p.add_argument('--min_diff', type=float, default=0.02,
                   help='Minimum difference between the two peaks to be '
                        'considered significant.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_fixel_mtr])
    assert_outputs_exist(parser, args, [args.out_differences,
                                        args.out_diff_mask,
                                        args.out_no_diff_mask,
                                        args.out_crossing_mask])

    peak_values_img = nib.load(args.in_fixel_mtr)
    peak_values = peak_values_img.get_fdata().astype(np.float32)

    mask = (peak_values[..., 0] >= args.min_mtr) & (peak_values[..., 1] >= args.min_mtr)

    peak_diff = np.where(mask, peak_values[..., 0] - peak_values[..., 1], 0)
    peak_diff = abs(peak_diff)

    diff_mask = peak_diff >= args.min_diff
    no_diff_mask = (peak_diff < args.min_diff) & (peak_values[..., 1] >= args.min_mtr)

    nib.save(nib.Nifti1Image(peak_diff, peak_values_img.affine),
             args.out_differences)
    nib.save(nib.Nifti1Image(diff_mask.astype(np.uint8),
                             peak_values_img.affine),
             args.out_diff_mask)
    nib.save(nib.Nifti1Image(no_diff_mask.astype(np.uint8),
                             peak_values_img.affine),
             args.out_no_diff_mask)
    nib.save(nib.Nifti1Image(mask.astype(np.uint8),
                             peak_values_img.affine),
             args.out_crossing_mask)


if __name__ == "__main__":
    main()
