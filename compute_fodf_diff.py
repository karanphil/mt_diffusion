#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import argparse
import logging

import nibabel as nib
import numpy as np

from dipy.data import get_sphere
from dipy.reconst.shm import sh_to_sf_matrix, order_from_ncoef

from scilpy.io.image import get_data_as_mask
from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('in_fodf_mt_off',
                   help='The DW image file to split.')

    p.add_argument('in_fodf_mt_on',
                   help='The DW image file to split.')

    p.add_argument('out_fodf')

    p.add_argument('--mask')

    p.add_argument('--sphere', default='repulsion724',
                   choices=['symmetric362', 'symmetric642', 'symmetric724',
                            'repulsion724', 'repulsion100', 'repulsion200'])

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_fodf_mt_off, args.in_fodf_mt_on])

    img_mt_off = nib.load(args.in_fodf_mt_off)
    img_mt_on = nib.load(args.in_fodf_mt_on)

    sh_mt_off = img_mt_off.get_fdata()
    sh_mt_on = img_mt_on.get_fdata()

    sphere = get_sphere(args.sphere)

    sh_order = order_from_ncoef(sh_mt_off.shape[-1])
    B, invB = sh_to_sf_matrix(sphere, sh_order)

    fodf_mt_off = np.dot(sh_mt_off, B)
    fodf_mt_off = fodf_mt_off.clip(min=0)

    fodf_mt_on = np.dot(sh_mt_on, B)
    fodf_mt_on = fodf_mt_on.clip(min=0)

    fodf_diff = np.where(fodf_mt_off >= 0, (fodf_mt_off - fodf_mt_on), 0)
    # fodf_diff = np.where(fodf_mt_off > 0, (fodf_mt_off - fodf_mt_on) / fodf_mt_off, 0)

    if args.mask:
        mask = get_data_as_mask(nib.load(args.mask), dtype=bool)
        reshaped_mask = np.repeat(mask, fodf_diff.shape[-1]).reshape(fodf_diff.shape)
        fodf_diff = np.where(np.squeeze(reshaped_mask), np.squeeze(fodf_diff), 0)

    sh_diff = np.dot(fodf_diff, invB)

    nib.save(nib.Nifti1Image(sh_diff, img_mt_off.affine), args.out_fodf)


if __name__ == "__main__":
    main()
