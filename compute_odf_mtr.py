#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import argparse
import logging

import nibabel as nib
import numpy as np

from dipy.data import get_sphere
from dipy.reconst.shm import sh_to_sf_matrix, order_from_ncoef, sf_to_sh

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

    p.add_argument('in_peaks_mt_off')

    p.add_argument('in_peaks_mt_on')

    p.add_argument('in_fixel_density_voxel_norm')

    p.add_argument('in_fixel_density_none_norm')

    p.add_argument('out_fodf')
    
    p.add_argument('out_peak_values')

    p.add_argument('out_peaks')

    p.add_argument('--mask')

    p.add_argument('--min_angle', default=10, type=float)

    p.add_argument('--rel_thr', default=0.1, type=float)

    p.add_argument('--abs_thr', default=1.5, type=float)

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
    assert_outputs_exist(parser, args, [args.out_fodf, args.out_peak_values,
                                        args.out_peaks])

    img_mt_off = nib.load(args.in_fodf_mt_off)
    img_mt_on = nib.load(args.in_fodf_mt_on)

    sh_mt_off = img_mt_off.get_fdata().astype(np.float32)
    sh_mt_on = img_mt_on.get_fdata().astype(np.float32)

    vol = nib.load(args.in_fixel_density_voxel_norm)
    fd_voxel = vol.get_fdata()
    fd_voxel = np.sum(fd_voxel, axis=-1)

    vol = nib.load(args.in_fixel_density_none_norm)
    fd_none = vol.get_fdata()
    fd_none = np.sum(fd_none, axis=-1)

    img_peaks_mt_off = nib.load(args.in_peaks_mt_off)
    peaks_mt_off = img_peaks_mt_off.get_fdata()

    img_peaks_mt_on = nib.load(args.in_peaks_mt_on)
    peaks_mt_on = img_peaks_mt_on.get_fdata()

    sphere = get_sphere(args.sphere)
    # subdivise?

    sh_order = order_from_ncoef(sh_mt_off.shape[-1])
    B, invB = sh_to_sf_matrix(sphere, sh_order)

    fodf_mt_off = np.dot(sh_mt_off, B)
    fodf_mt_off = fodf_mt_off.clip(min=0)

    fodf_mt_on = np.dot(sh_mt_on, B)
    fodf_mt_on = fodf_mt_on.clip(min=0)

    nb_peaks = 5
    nb_vertices = len(sphere.vertices)

    mtr_sphere = np.zeros(peaks_mt_off.shape[0:3] + (nb_vertices,))
    mtr_peaks = np.zeros(peaks_mt_off.shape[0:3] + (nb_peaks,))
    new_peaks = np.zeros(peaks_mt_off.shape)
    missed = 0
    for i in range(mtr_sphere.shape[0]):
        for j in range(mtr_sphere.shape[1]):
            for k in range(mtr_sphere.shape[2]):
                peak_indices = [0, 1, 2, 3, 4]
                for l in range(nb_peaks):
                    if fd_voxel[i,j,k,l] > args.rel_thr and fd_none[i,j,k,l] > args.abs_thr:
                        found_peak = False
                        for m in peak_indices:
                            flip = 1
                            # Compute angle between peaks
                            angle_between = np.rad2deg(np.arccos(np.clip(np.dot(peaks_mt_off[i,j,k, 3*l:3*(l+1)],peaks_mt_on[i,j,k, 3*m:3*(m+1)]), -1.0, 1.0)))
                            angle_between_flip = np.rad2deg(np.arccos(np.clip(np.dot(peaks_mt_off[i,j,k, 3*l:3*(l+1)],-1*peaks_mt_on[i,j,k, 3*m:3*(m+1)]), -1.0, 1.0)))
                            if angle_between_flip < angle_between:
                                angle_between = angle_between_flip
                                flip = -1
                            if angle_between <= args.min_angle:
                                peak_indices.remove(m)
                                found_peak = True
                                mean_peak = (peaks_mt_off[i,j,k, 3*l:3*(l+1)] + flip * peaks_mt_on[i,j,k, 3*m:3*(m+1)]) / 2
                                new_peaks[i,j,k, 3*l:3*(l+1)] = mean_peak / np.linalg.norm(mean_peak)
                                vector_mt_off = sphere.find_closest(peaks_mt_off[i,j,k, 3*l:3*(l+1)])
                                vector_mt_on = sphere.find_closest(peaks_mt_on[i,j,k, 3*m:3*(m+1)])
                                vector = sphere.find_closest(mean_peak)
                                mtr_sphere[i,j,k,vector] = (fodf_mt_off[i,j,k,vector_mt_off] - fodf_mt_on[i,j,k,vector_mt_on]) / fodf_mt_off[i,j,k,vector_mt_off]
                                mtr_peaks[i,j,k,l] = (fodf_mt_off[i,j,k,vector_mt_off] - fodf_mt_on[i,j,k,vector_mt_on]) / fodf_mt_off[i,j,k,vector_mt_off]
                        if not found_peak:
                            missed += 1
                            logging.warning(f'No matching peak found for voxel {i},{j},{k} peak {l}')
    logging.warning(f'Total missed peaks: {missed}')
    mtr_sphere = np.clip(mtr_sphere, a_min=0, a_max=None)
    mtr_peaks = np.clip(mtr_peaks, a_min=0, a_max=None)

    if args.mask:
        mask = get_data_as_mask(nib.load(args.mask), dtype=bool)
        reshaped_mask = np.repeat(mask, mtr_sphere.shape[-1]).reshape(mtr_sphere.shape)
        mtr_sphere = np.where(np.squeeze(reshaped_mask), np.squeeze(mtr_sphere), 0)
        reshaped_mask = np.repeat(mask, mtr_peaks.shape[-1]).reshape(mtr_peaks.shape)
        mtr_peaks = np.where(np.squeeze(reshaped_mask), np.squeeze(mtr_peaks), 0)

    mtr_sh = sf_to_sh(mtr_sphere, sphere, sh_order_max=6, smooth=0.1)

    nib.save(nib.Nifti1Image(mtr_sh, img_mt_off.affine), args.out_fodf)
    nib.save(nib.Nifti1Image(mtr_peaks, img_mt_off.affine), args.out_peak_values)
    nib.save(nib.Nifti1Image(new_peaks.astype(np.float32), img_mt_off.affine), args.out_peaks)

    # fodf_diff = np.where(fodf_mt_off >= 0, (fodf_mt_off - fodf_mt_on), 0)
    # # fodf_diff = np.where(fodf_mt_off > 0, (fodf_mt_off - fodf_mt_on) / fodf_mt_off, 0)

    # # fodf_ratio = np.where(fodf_mt_off > args.thr, (fodf_diff / fodf_mt_off), 0)
    # # fodf_ratio = np.where(fodf_mt_off * 0.2 > args.thr, (fodf_diff / fodf_mt_off), 0)
    # fodf_ratio = np.where(fodf_mt_off * 0.2 > args.thr, (fodf_diff / fodf_mt_off), 0)

    # if args.mask:
    #     mask = get_data_as_mask(nib.load(args.mask), dtype=bool)
    #     reshaped_mask = np.repeat(mask, fodf_ratio.shape[-1]).reshape(fodf_ratio.shape)
    #     fodf_ratio = np.where(np.squeeze(reshaped_mask), np.squeeze(fodf_ratio), 0)

    # sh_ratio = np.dot(fodf_ratio, invB)

    # nib.save(nib.Nifti1Image(sh_ratio, img_mt_off.affine), args.out_fodf2)


if __name__ == "__main__":
    main()
