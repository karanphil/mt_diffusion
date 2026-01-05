#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to compute the average fixel-MTR for each voxel. Averages the MTR of
each peak (fixel) in a voxel, weighted by their peak fraction. The peak
fractions are computed from the peak values provided from the MT-off fODFs.
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

    p.add_argument('in_mtr_peak_values',
                   help='MTR peak values file, obtained from '
                        'compute_fixel_mtr.py.')

    p.add_argument('in_peak_values',
                   help='Peak values files, obtained from scil_fodf_metrics. '
                        'Should come from the MT-off fODFs. Is used to '
                        'compute the peak fractions.')

    p.add_argument('out_average_fixel_mtr',
                   help='Output average fixel-MTR file.')
    
    p.add_argument('--out_peak_fractions',
                   help='If provided, will output the computed peak fractions '
                        'file.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_mtr_peak_values,
                                 args.in_peak_values])
    assert_outputs_exist(parser, args, [args.out_average_fixel_mtr])

    mtr_peak_values_img = nib.load(args.in_mtr_peak_values)
    mtr_peak_values = mtr_peak_values_img.get_fdata().astype(np.float32)

    peak_values_img = nib.load(args.in_peak_values)
    peak_values = peak_values_img.get_fdata().astype(np.float32)

    peak_norm = np.sum(peak_values, axis=-1)

    for i in range(peak_values.shape[-1]):
        peak_values[..., i] =  np.where(peak_norm > 0,
                                        peak_values[..., i] / peak_norm, 0)

    if args.out_peak_fractions:
        nib.save(nib.Nifti1Image(peak_values, peak_values_img.affine),
                 args.out_peak_fractions)

    peak_values[mtr_peak_values == 0] = 0
    peak_norm = np.sum(peak_values, axis=-1)
    peak_values[peak_norm == 0] = np.nan
    mtr = np.average(mtr_peak_values, weights=peak_values, axis=-1)

    nib.save(nib.Nifti1Image(mtr, peak_values_img.affine),
             args.out_average_fixel_mtr)


if __name__ == "__main__":
    main()
