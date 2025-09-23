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

    p.add_argument('in_images', nargs='+',
                   help='Input images to average.')

    p.add_argument('out_average')

    p.add_argument('--median', action='store_true',
                   help='Compute the median instead of the mean.')
    
    p.add_argument('--exclude_zeros', action='store_true',
                   help='Exclude zeros from the average computation.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, args.in_images)
    assert_outputs_exist(parser, args, [args.out_average])

    for idx, in_file in enumerate(args.in_images):
        img = nib.load(in_file)
        data = img.get_fdata(dtype=np.float32)

        if idx == 0:
            all_data = np.zeros(data.shape + (len(args.in_images),), dtype=np.float32)

        all_data[..., idx] = data

    if args.median:
        out_data = np.median(all_data, axis=-1)
    else:
        if args.exclude_zeros:
            # This option works weird, it produces noise on the sides of bundles where few subjects have data
            out_data = np.divide(np.sum(all_data, axis=-1),
                                 np.count_nonzero(all_data, axis=-1),
                                 out=np.zeros_like(all_data[..., 0]),
                                 where=np.count_nonzero(all_data, axis=-1) != 0)
        else:
            out_data = np.mean(all_data, axis=-1)

    nib.save(nib.Nifti1Image(out_data, img.affine), args.out_average)


if __name__ == "__main__":
    main()
