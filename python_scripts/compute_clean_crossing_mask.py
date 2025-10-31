#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to compute a clean crossing mask based on the number of fibers,
number of bundles and fixel density maps.
"""

import argparse
import logging

import nibabel as nib
import numpy as np

from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('in_nufo',
                   help='Number of fiber orientations per voxel. This should '
                        'be the result of scil_fodf_metrics, named as '
                        'nufo.nii.gz')

    p.add_argument('in_nb_bundles',
                   help='Number of bundles per voxel. This should be the '
                        'result of scil_bundle_fixel_analysis, named as '
                        'nb_bundles_per_voxel_voxel-norm.nii.gz')

    p.add_argument('in_fixel_density',
                   help='Fixel density map for each bundle. This should be '
                        'the result of scil_bundle_fixel_analysis, named '
                        'as fixel_density_maps_voxel-norm.nii.gz')

    p.add_argument('out_mask',
                   help='Output clean crossing mask.')

    p.add_argument('--thr', default=0.7, type=float,
                   help='Threshold on the fixel density map to consider a '
                        'voxel as a clean crossing. Default is 0.7.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_nufo, args.in_nb_bundles,
                                 args.in_fixel_density])

    vol = nib.load(args.in_nufo)
    nufo = vol.get_fdata()
    vol = nib.load(args.in_nb_bundles)
    nb_bundles = vol.get_fdata()
    vol = nib.load(args.in_fixel_density)
    fd = vol.get_fdata()

    # If the fixel density of any bundle is above the threshold (for the 
    # first peak), we don't keep it. This ensures crossing with maximum 70-30
    # peak fraction.
    fd_mask = np.all(fd[:, :, :, 0, :] <= args.thr, axis=-1)

    mask = (nufo == 2) & (nb_bundles == 2) & fd_mask

    nib.save(nib.Nifti1Image(mask.astype(np.uint8), vol.affine), args.out_mask)


if __name__ == "__main__":
    main()
