import argparse
import nibabel as nib
import numpy as np
from pathlib import Path

from modules.utils import (compute_single_fiber_averages,
                           fit_single_fiber_results,
                           correct_measure)
from modules.io import save_txt, plot_means, save_masks_by_angle_bins

from scilpy.io.utils import (add_overwrite_arg, add_processes_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('in_peaks',
                   help='Path of the fODF peaks.')
    p.add_argument('in_fa',
                   help='Path of the FA.')
    p.add_argument('in_wm_mask',
                   help='Path of the WM mask.')
    p.add_argument('out_folder',
                   help='Path of the output folder for txt, png, masks and measures.')
    
    p.add_argument('--in_e1',
                   help='Path to the principal eigenvector of DTI.')
    p.add_argument('--in_nufo',
                   help='Path to the NuFO.')
    p.add_argument('--in_mtr',
                   help='Path to the MTR.')
    p.add_argument('--in_ihmtr',
                   help='Path to the ihMTR.')
    p.add_argument('--in_mtsat',
                   help='Path to the MTsat.')
    p.add_argument('--in_ihmtsat',
                   help='Path to the ihMTsat.')

    p.add_argument('--files_basename',
                   help='Basename of all the saved txt or png files.')

    p.add_argument('--fa_thr', default=0.5,
                   help='Value of FA threshold [%(default)s].')
    p.add_argument('--bin_width', default=10,
                   help='Value of the bin width.')
    p.add_argument('--poly_order', default=8,
                   help='Order of the polynome to fit [%(default)s].')
    add_overwrite_arg(p)
    # add_processes_arg(p)
    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.files_basename:
        files_basename = args.files_basename
    else:
        files_basename = "_" + str(args.fa_thr) + "_fa_thr_" \
            + str(args.bin_width) + "_bin_width"
        
    out_folder = Path(args.out_folder)

    # Load the data
    peaks_img = nib.load(args.in_peaks)
    fa_img = nib.load(args.in_fa)
    wm_mask_img = nib.load(args.in_wm_mask)

    peaks = peaks_img.get_fdata()
    fa = fa_img.get_fdata()
    wm_mask = wm_mask_img.get_fdata()

    affine = peaks_img.affine

    if args.in_e1:
        e1_img = nib.load(args.in_e1)
        e1 = e1_img.get_fdata()
    else:
        e1 = peaks
    if args.in_nufo:
        nufo_img = nib.load(args.in_nufo)
        nufo = nufo_img.get_fdata()
    else:
        nufo = None
    if args.in_mtr:
        mtr_img = nib.load(args.in_mtr)
        mtr = mtr_img.get_fdata()
    else:
        mtr = None
    if args.in_ihmtr:
        ihmtr_img = nib.load(args.in_ihmtr)
        ihmtr = ihmtr_img.get_fdata()
    else:
        ihmtr = None
    if args.in_mtsat:
        mtsat_img = nib.load(args.in_mtsat)
        mtsat = mtsat_img.get_fdata()
    else:
        mtsat = None
    if args.in_ihmtsat:
        ihmtsat_img = nib.load(args.in_ihmtsat)
        ihmtsat = ihmtsat_img.get_fdata()
    else:
        ihmtsat = None

    print("Computing single fiber averages.")
    bins, mtr_means, ihmtr_means, mtsat_means, ihmtsat_means, \
        nb_voxels = compute_single_fiber_averages(e1, fa, wm_mask, affine,
                                                  mtr=mtr, ihmtr=ihmtr,
                                                  mtsat=mtsat,
                                                  ihmtsat=ihmtsat,
                                                  nufo=nufo,
                                                  bin_width=args.bin_width,
                                                  fa_thr=args.fa_thr)

    print("Saving results as txt files.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "single_fiber_mtr_ihmtr_results" + files_basename + ".txt"
        output_path = out_folder / output_name
        save_txt(bins, mtr_means, ihmtr_means, nb_voxels,
                 str(output_path), input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "single_fiber_mtsat_ihmtsat_results" + files_basename + ".txt"
        output_path = out_folder / output_name
        save_txt(bins, mtsat_means, ihmtsat_means, nb_voxels,
                 str(output_path), input_dtype="sats")
        
    print("Fitting the results.")
    if args.in_mtr and args.in_ihmtr:
        mtr_poly = fit_single_fiber_results(bins, mtr_means,
                                            poly_order=args.poly_order)
        ihmtr_poly = fit_single_fiber_results(bins, ihmtr_means,
                                              poly_order=args.poly_order)
    if args.in_mtsat and args.in_ihmtsat:
        mtsat_poly = fit_single_fiber_results(bins, mtsat_means,
                                            poly_order=args.poly_order)
        ihmtsat_poly = fit_single_fiber_results(bins, ihmtsat_means,
                                              poly_order=args.poly_order)

    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "single_fiber_mtr_ihmtr_plot" + files_basename + ".png"
        output_path = out_folder / output_name
        plot_means(bins, mtr_means, ihmtr_means, nb_voxels, str(output_path),
                   mtr_poly, ihmtr_poly, input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "single_fiber_mtsat_ihmtsat_plot" + files_basename + ".png"
        output_path = out_folder / output_name
        plot_means(bins, mtsat_means, ihmtsat_means, nb_voxels,
                   str(output_path), mtsat_poly, ihmtsat_poly,
                   input_dtype="sats")

    print("Saving single fiber masks.")
    save_masks_by_angle_bins(e1, fa, wm_mask, affine,
                             out_folder, nufo=nufo,
                             bin_width=args.bin_width, fa_thr=args.fa_thr)
    
    print("Correcting measures.")
    if args.in_mtr:
        corrected_mtr = correct_measure(peaks, mtr, affine, wm_mask, mtr_poly,
                                        peak_frac_thr=0.1)
        corrected_name = "mtr_corrected.nii.gz"
        corrected_path = out_folder / corrected_name
        nib.save(nib.Nifti1Image(corrected_mtr, affine), corrected_path)
    if args.in_ihmtr:
        corrected_ihmtr = correct_measure(peaks, ihmtr, affine, wm_mask,
                                          ihmtr_poly, peak_frac_thr=0.1)
        corrected_name = "ihmtr_corrected.nii.gz"
        corrected_path = out_folder / corrected_name
        nib.save(nib.Nifti1Image(corrected_ihmtr, affine), corrected_path)
    if args.in_mtsat:
        corrected_mtsat = correct_measure(peaks, mtsat, affine, wm_mask,
                                          mtsat_poly, peak_frac_thr=0.1)
        corrected_name = "mtsat_corrected.nii.gz"
        corrected_path = out_folder / corrected_name
        nib.save(nib.Nifti1Image(corrected_mtsat, affine), corrected_path)
    if args.in_ihmtr:
        corrected_ihmtsat = correct_measure(peaks, ihmtsat, affine, wm_mask,
                                          ihmtsat_poly, peak_frac_thr=0.1)
        corrected_name = "ihmtsat_corrected.nii.gz"
        corrected_path = out_folder / corrected_name
        nib.save(nib.Nifti1Image(corrected_ihmtsat, affine), corrected_path)


if __name__ == "__main__":
    main()
