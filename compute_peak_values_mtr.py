#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import argparse
import logging

import nibabel as nib
import numpy as np

from dipy.data import get_sphere
from dipy.reconst.shm import sf_to_sh

from scilpy.io.image import get_data_as_mask
from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('in_peak_values_mt_off')

    p.add_argument('in_peak_values_mt_on')

    p.add_argument('in_peaks')

    p.add_argument('out_peak_values')

    p.add_argument('out_fodf')

    p.add_argument('--mask')

    # p.add_argument('--thr', default=0, type=float)

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

    assert_inputs_exist(parser, [args.in_peak_values_mt_off, args.in_peak_values_mt_on, args.in_peaks])

    img_mt_off = nib.load(args.in_peak_values_mt_off)
    img_mt_on = nib.load(args.in_peak_values_mt_on)

    pvalues_mt_off = img_mt_off.get_fdata()
    pvalues_mt_on = img_mt_on.get_fdata()

    peak_values_diff = np.where(pvalues_mt_off > 0, pvalues_mt_off - pvalues_mt_on, 0)
    peak_values_mtr = np.where(pvalues_mt_off > 0, peak_values_diff / pvalues_mt_off, 0)

    if args.mask:
        mask = get_data_as_mask(nib.load(args.mask), dtype=bool)
        reshaped_mask = np.repeat(mask, peak_values_mtr.shape[-1]).reshape(peak_values_mtr.shape)
        peak_values_mtr = np.where(np.squeeze(reshaped_mask), np.squeeze(peak_values_mtr), 0)

    nib.save(nib.Nifti1Image(peak_values_mtr, img_mt_off.affine), args.out_peak_values)

    sphere = get_sphere(args.sphere)
    sphere = sphere.subdivide(1)

    nb_peaks = 3
    nb_vertices = len(sphere.vertices)

    img_peaks = nib.load(args.in_peaks)
    peaks = img_peaks.get_fdata()

    mtr_sphere = np.zeros(peak_values_mtr.shape[0:3] + (nb_vertices,))
    for i in range(mtr_sphere.shape[0]):
        for j in range(mtr_sphere.shape[1]):
            for k in range(mtr_sphere.shape[2]):
                for l in range(nb_peaks):
                    if peak_values_mtr[i,j,k,l] != 0:
                        vector = sphere.find_closest(peaks[i,j,k, 3*l:3*(l+1)])
                        mtr_sphere[i,j,k,vector] = peak_values_mtr[i,j,k,l]
    
    mtr_sh = sf_to_sh(mtr_sphere, sphere, sh_order_max=8, smooth=0.1)

    nib.save(nib.Nifti1Image(mtr_sh, img_mt_off.affine), args.out_fodf)


if __name__ == "__main__":
    main()
