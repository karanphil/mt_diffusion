#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to compute the fixel MTR from MT-on and MT-off fODF amplitudes.
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
                   help='The input fODF MT-off file.')

    p.add_argument('in_fodf_mt_on',
                   help='The input fODF MT-on file.')

    p.add_argument('in_peaks_mt_off',
                   help='The input peaks MT-off file.')

    p.add_argument('in_peaks_mt_on',
                   help='The input peaks MT-on file.')

    p.add_argument('in_fixel_density_voxel_norm',
                   help='Fixel density map normalized by voxel for each '
                        'bundle. This should be the result of '
                        'scil_bundle_fixel_analysis, named '
                        'as fixel_density_maps_voxel-norm.nii.gz')

    p.add_argument('in_fixel_density_none_norm',
                   help='Fixel density map not normalized for each '
                        'bundle. This should be the result of '
                        'scil_bundle_fixel_analysis, named '
                        'as fixel_density_maps_none-norm.nii.gz')

    p.add_argument('out_fodf',
                   help='Output FODF MTR. IMPORTANT: The amplitudes of the '
                        'fODFs have no meaning.')
    
    p.add_argument('out_peak_values',
                   help='Output fixel-MTR.')

    p.add_argument('out_peaks',
                   help='Output peaks with MT-on and MT-off averaged.')

    p.add_argument('--mask',
                   help='Optional mask to limit the computation.')

    p.add_argument('--min_angle', default=10, type=float,
                   help='Minimum angle (in degrees) to consider two peaks '
                        'as being the same direction. Default is 10 degrees.')

    p.add_argument('--rel_thr', default=0.01, type=float,
                   help='Relative threshold on the voxel-normalized fixel '
                        'density map to consider a fixel for MTR '
                        'computation. Default is 0.01 (1%% of max).')

    p.add_argument('--abs_thr', default=0.0, type=float,
                   help='Absolute threshold on the non-normalized fixel '
                        'density map to consider a fixel for MTR '
                        'computation. Default is 0.')

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

    if args.mask:
        mask = get_data_as_mask(nib.load(args.mask), dtype=bool)
    else:
        mask = np.ones(sh_mt_off.shape[0:3], dtype=bool)

    sphere = get_sphere(args.sphere)

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
                    if fd_voxel[i,j,k,l] > args.rel_thr and fd_none[i,j,k,l] > args.abs_thr and mask[i,j,k]:
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
                            logging.debug(f'No matching peak found for voxel {i},{j},{k} peak {l}')
    logging.warning(f'Total missed peaks: {missed}')
    mtr_sphere = np.clip(mtr_sphere, a_min=0, a_max=1)
    mtr_peaks = np.clip(mtr_peaks, a_min=0, a_max=1)

    mtr_sh = sf_to_sh(mtr_sphere * 50, sphere, sh_order_max=6, smooth=0.1)

    nib.save(nib.Nifti1Image(mtr_sh, img_mt_off.affine), args.out_fodf)
    nib.save(nib.Nifti1Image(mtr_peaks, img_mt_off.affine), args.out_peak_values)
    nib.save(nib.Nifti1Image(new_peaks.astype(np.float32), img_mt_off.affine), args.out_peaks)


if __name__ == "__main__":
    main()
