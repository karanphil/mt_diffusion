#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to remove NaN values from a DWI NIfTI file by replacing them with zeros.
"""

import argparse
import logging

import nibabel as nib
import numpy as np

from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('in_dwi', help='Input DWI file (nifti).')

    p.add_argument('out_dwi', help='Output DWI file (nifti).')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_dwi])
    assert_outputs_exist(parser, args, [args.out_dwi])

    vol = nib.load(args.in_dwi)
    dwi = vol.get_fdata().astype(np.float32)
    new_dwi = np.where(np.isnan(dwi), 0, dwi)
    nib.save(nib.Nifti1Image(new_dwi.astype(np.float32), vol.affine), args.out_dwi)


if __name__ == "__main__":
    main()