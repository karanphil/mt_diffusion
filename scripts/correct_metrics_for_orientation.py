import argparse
import nibabel as nib
import numpy as np
from pathlib import Path

from modules.utils import (compute_single_fiber_averages,
                           compute_crossing_fibers_averages,
                           fit_single_fiber_results,
                           correct_measure)
from modules.io import save_txt, plot_means, save_masks_by_angle_bins

from scilpy.io.utils import (add_overwrite_arg)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('in_peaks',
                   help='Path of the fODF peaks.')
    p.add_argument('in_peak_values',
                   help='Path of the fODF peak values.')
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
    p.add_argument('--bin_width', default=1,
                   help='Value of the bin width [%(default)s].')
    p.add_argument('--norm_thr', default=0.7,
                   help='Value of the norm threshold for selecting 2 fibers [%(default)s].')
    p.add_argument('--poly_order', default=10,
                   help='Order of the polynome to fit [%(default)s].')
    add_overwrite_arg(p)
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
    peak_values_img = nib.load(args.in_peak_values)
    fa_img = nib.load(args.in_fa)
    wm_mask_img = nib.load(args.in_wm_mask)

    peaks = peaks_img.get_fdata()
    peak_values = peak_values_img.get_fdata()
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
                   mt_poly=mtr_poly, ihmt_poly=ihmtr_poly,
                   input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "single_fiber_mtsat_ihmtsat_plot" + files_basename + ".png"
        output_path = out_folder / output_name
        plot_means(bins, mtsat_means, ihmtsat_means, nb_voxels,
                   str(output_path), mt_poly=mtsat_poly,
                   ihmt_poly=ihmtsat_poly, input_dtype="sats")

    print("Saving single fiber masks.")
    save_masks_by_angle_bins(e1, fa, wm_mask, affine,
                             out_folder, nufo=nufo,
                             bin_width=10, fa_thr=args.fa_thr)
    
    # print("Computing crossing fibers average.")
    # bins_2f, mtr_2f_means, ihmtr_2f_means, mtsat_2f_means, ihmtsat_2f_means, \
    #     nb_voxels_2f = compute_crossing_fibers_averages(peaks, wm_mask, affine,
    #                                                     nufo, mtr=mtr,
    #                                                     ihmtr=ihmtr,
    #                                                     mtsat=mtsat,
    #                                                     ihmtsat=ihmtsat,
    #                                                     bin_width=10,
    #                                                     norm_thr=args.norm_thr)
    
    # mtr_2f_means_diag = np.diagonal(mtr_2f_means)
    # ihmtr_2f_means_diag = np.diagonal(ihmtr_2f_means)
    # mtsat_2f_means_diag = np.diagonal(mtsat_2f_means)
    # ihmtsat_2f_means_diag = np.diagonal(ihmtsat_2f_means)
    # nb_voxels_2f_diag = np.diagonal(nb_voxels_2f)
    
    # print("Saving results as plots.")
    # if args.in_mtr and args.in_ihmtr:
    #     output_name = "double_fibers_mtr_ihmtr_diagonal_plot" + files_basename + ".png"
    #     output_path = out_folder / output_name
    #     plot_means(bins_2f, mtr_2f_means_diag, ihmtr_2f_means_diag,
    #                nb_voxels_2f_diag, str(output_path), input_dtype="ratios")
    # if args.in_mtsat and args.in_ihmtsat:
    #     output_name = "double_fibers_mtsat_ihmtsat_diagonal_plot" + files_basename + ".png"
    #     output_path = out_folder / output_name
    #     plot_means(bins_2f, mtsat_2f_means_diag, ihmtsat_2f_means_diag,
    #                nb_voxels_2f_diag, str(output_path), input_dtype="sats")
    
    print("Correcting measures.")
    if args.in_mtr:
        corrected_mtr = correct_measure(peaks, peak_values, mtr, affine,
                                        wm_mask, mtr_poly,
                                        peak_frac_thr=0.0)
        corrected_name = "mtr_corrected.nii.gz"
        corrected_path = out_folder / corrected_name
        nib.save(nib.Nifti1Image(corrected_mtr, affine), corrected_path)
    if args.in_ihmtr:
        corrected_ihmtr = correct_measure(peaks, peak_values, ihmtr, affine,
                                          wm_mask, ihmtr_poly,
                                          peak_frac_thr=0.0)
        corrected_name = "ihmtr_corrected.nii.gz"
        corrected_path = out_folder / corrected_name
        nib.save(nib.Nifti1Image(corrected_ihmtr, affine), corrected_path)
    if args.in_mtsat:
        corrected_mtsat = correct_measure(peaks, peak_values, mtsat, affine,
                                          wm_mask, mtsat_poly,
                                          peak_frac_thr=0.0)
        corrected_name = "mtsat_corrected.nii.gz"
        corrected_path = out_folder / corrected_name
        nib.save(nib.Nifti1Image(corrected_mtsat, affine), corrected_path)
    if args.in_ihmtr:
        corrected_ihmtsat = correct_measure(peaks, peak_values, ihmtsat,
                                            affine, wm_mask, ihmtsat_poly,
                                            peak_frac_thr=0.0)
        corrected_name = "ihmtsat_corrected.nii.gz"
        corrected_path = out_folder / corrected_name
        nib.save(nib.Nifti1Image(corrected_ihmtsat, affine), corrected_path)

    print("Computing single fiber averages on corrected data.")
    bins, mtr_cr_means, ihmtr_cr_means, mtsat_cr_means, ihmtsat_cr_means, \
        nb_voxels = compute_single_fiber_averages(e1, fa, wm_mask, affine,
                                                  mtr=corrected_mtr,
                                                  ihmtr=corrected_ihmtr,
                                                  mtsat=corrected_mtsat,
                                                  ihmtsat=corrected_ihmtsat,
                                                  nufo=nufo,
                                                  bin_width=args.bin_width,
                                                  fa_thr=args.fa_thr)

    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "single_fiber_corrected_mtr_ihmtr_plot" + files_basename + ".png"
        output_path = out_folder / output_name
        plot_means(bins, mtr_means, ihmtr_means, nb_voxels,
                   str(output_path), mt_cr_means=mtr_cr_means,
                   ihmt_cr_means=ihmtr_cr_means, mt_poly=mtr_poly,
                   ihmt_poly=ihmtr_poly, input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "single_fiber_corrected_mtsat_ihmtsat_plot" + files_basename + ".png"
        output_path = out_folder / output_name
        plot_means(bins, mtsat_means, ihmtsat_means, nb_voxels,
                   str(output_path), mt_cr_means=mtsat_cr_means,
                   ihmt_cr_means=ihmtsat_cr_means, mt_poly=mtsat_poly,
                   ihmt_poly=ihmtsat_poly, input_dtype="sats")
        
    # print("Computing crossing fibers average on corrected data.")
    # bins_2f, mtr_cr_2f_means, ihmtr_cr_2f_means, mtsat_cr_2f_means, ihmtsat_cr_2f_means, \
    #     nb_voxels_2f = compute_crossing_fibers_averages(peaks, wm_mask, affine,
    #                                                     nufo, mtr=corrected_mtr,
    #                                                     ihmtr=corrected_ihmtr,
    #                                                     mtsat=corrected_mtsat,
    #                                                     ihmtsat=corrected_ihmtsat,
    #                                                     bin_width=10,
    #                                                     norm_thr=args.norm_thr)
    
    # mtr_cr_2f_means_diag = np.diagonal(mtr_cr_2f_means)
    # ihmtr_cr_2f_means_diag = np.diagonal(ihmtr_cr_2f_means)
    # mtsat_cr_2f_means_diag = np.diagonal(mtsat_cr_2f_means)
    # ihmtsat_cr_2f_means_diag = np.diagonal(ihmtsat_cr_2f_means)
    # nb_voxels_2f_diag = np.diagonal(nb_voxels_2f)
    
    # print("Saving results as plots.")
    # if args.in_mtr and args.in_ihmtr:
    #     output_name = "double_fibers_corrected_mtr_ihmtr_diagonal_plot" + files_basename + ".png"
    #     output_path = out_folder / output_name
    #     plot_means(bins_2f, mtr_2f_means_diag, ihmtr_2f_means_diag,
    #                nb_voxels_2f_diag, str(output_path), 
    #                mt_cr_means=mtr_cr_2f_means_diag,
    #                ihmt_cr_means=ihmtr_cr_2f_means_diag, input_dtype="ratios")
    # if args.in_mtsat and args.in_ihmtsat:
    #     output_name = "double_fibers_corrected_mtsat_ihmtsat_diagonal_plot" + files_basename + ".png"
    #     output_path = out_folder / output_name
    #     plot_means(bins_2f, mtsat_2f_means_diag, ihmtsat_2f_means_diag,
    #                nb_voxels_2f_diag, str(output_path),
    #                mt_cr_means=mtsat_cr_2f_means_diag,
    #                ihmt_cr_means=ihmtsat_cr_2f_means_diag, input_dtype="sats")

if __name__ == "__main__":
    main()
