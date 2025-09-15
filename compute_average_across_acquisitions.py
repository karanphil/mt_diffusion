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
            avg_data = np.zeros(data.shape, dtype=np.float32)

        avg_data += data 
        
    avg_data /= len(args.in_images)

    nib.save(nib.Nifti1Image(avg_data, img.affine), args.out_average)


if __name__ == "__main__":
    main()
