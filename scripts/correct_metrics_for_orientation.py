import argparse
import nibabel as nib
import numpy as np
from pathlib import Path

from modules.utils import (compute_single_fiber_averages,
                           compute_crossing_fibers_averages,
                           fit_single_fiber_results,
                           correct_measure, compute_poor_ihmtr,
                           analyse_delta_m_max, analyse_3_crossing_fibers_averages)
from modules.io import (save_txt, plot_means, save_masks_by_angle_bins,
                        plot_multiple_means, plot_delta_m_max, plot_different_bins_means,
                        plot_3d_means)

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
    p.add_argument('in_nufo',
                   help='Path to the NuFO.')
    p.add_argument('out_folder',
                   help='Path of the output folder for txt, png, masks and measures.')
    
    p.add_argument('--in_e1',
                   help='Path to the principal eigenvector of DTI.')
    p.add_argument('--in_v1',
                   help='Path to the principal eigenvalue of DTI.')
    p.add_argument('--in_mtr',
                   help='Path to the MTR.')
    p.add_argument('--in_ihmtr',
                   help='Path to the ihMTR.')
    p.add_argument('--in_mtsat',
                   help='Path to the MTsat.')
    p.add_argument('--in_ihmtsat',
                   help='Path to the ihMTsat.')
    p.add_argument('--in_roi',
                   help='Path to the ROI for single fiber analysis.')
    p.add_argument('--in_bundles_folder',
                   help='Path to the bundles folder.')

    p.add_argument('--files_basename',
                   help='Basename of all the saved txt or png files.')

    p.add_argument('--fa_thr', default=0.5,
                   help='Value of FA threshold [%(default)s].')
    p.add_argument('--bin_width', default=1,
                   help='Value of the bin width for the whole brain [%(default)s].')
    p.add_argument('--bin_width_mask', default=3,
                   help='Value of the bin width inside the mask [%(default)s].')
    p.add_argument('--bin_width_bundles', default=5,
                   help='Value of the bin width inside bundles [%(default)s].')
    # p.add_argument('--frac_thr', default=0.4,
    #                help='Value of the fraction threshold for selecting 2 fibers [%(default)s].')
    p.add_argument('--min_frac_thr', default=0.1,
                   help='Value of the minimal fraction threshold for selecting peaks to correct [%(default)s].')
    p.add_argument('--min_nb_voxels', default=30, type=int,
                   help='Value of the minimal number of voxels per bin [%(default)s].')
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
    nufo_img = nib.load(args.in_nufo)


    peaks = peaks_img.get_fdata()
    peak_values = peak_values_img.get_fdata()
    fa = fa_img.get_fdata()
    wm_mask = wm_mask_img.get_fdata()
    nufo = nufo_img.get_fdata()

    affine = peaks_img.affine

    if args.in_e1:
        e1_img = nib.load(args.in_e1)
        e1 = e1_img.get_fdata()
    else:
        e1 = peaks
    if args.in_v1:
        v1_img = nib.load(args.in_v1)
        v1 = v1_img.get_fdata()
    else:
        v1 = peak_values
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
    if args.in_roi:
        roi_img = nib.load(args.in_roi)
        roi = roi_img.get_fdata()
    else:
        roi = None
    if args.in_bundles_folder:
        bundles_folder = Path(args.in_bundles_folder)
        bundles_list = list(bundles_folder.iterdir())
        nb_bundles = len(bundles_list)
        bundles = np.ndarray((wm_mask.shape) + (nb_bundles,))
        for i in range(nb_bundles):
            bundles[..., i] = nib.load(str(bundles_list[i])).get_fdata()
    else:
        bundles = None
        nb_bundles = 0

#----------------------------Whole brain----------------------------------------
    print("Computing single fiber averages for whole brain.")
    w_brain_results = compute_single_fiber_averages(e1, fa,
                                                    wm_mask,
                                                    affine,
                                                    mtr=mtr,
                                                    ihmtr=ihmtr,
                                                    mtsat=mtsat,
                                                    ihmtsat=ihmtsat,
                                                    nufo=nufo,
                                                    bin_width=args.bin_width,
                                                    fa_thr=args.fa_thr,
                                                    min_nb_voxels=args.min_nb_voxels)
    
    print("Computing single fiber averages for whole brain.")
    w_brain_results_3 = compute_single_fiber_averages(e1, fa,
                                                    wm_mask,
                                                    affine,
                                                    mtr=mtr,
                                                    ihmtr=ihmtr,
                                                    mtsat=mtsat,
                                                    ihmtsat=ihmtsat,
                                                    nufo=nufo,
                                                    bin_width=3,
                                                    fa_thr=args.fa_thr,
                                                    min_nb_voxels=args.min_nb_voxels)
    
    print("Computing single fiber averages for whole brain.")
    w_brain_results_5 = compute_single_fiber_averages(e1, fa,
                                                    wm_mask,
                                                    affine,
                                                    mtr=mtr,
                                                    ihmtr=ihmtr,
                                                    mtsat=mtsat,
                                                    ihmtsat=ihmtsat,
                                                    nufo=nufo,
                                                    bin_width=5,
                                                    fa_thr=args.fa_thr,
                                                    min_nb_voxels=args.min_nb_voxels)
    
    print("Computing single fiber averages for whole brain.")
    w_brain_results_10 = compute_single_fiber_averages(e1, fa,
                                                    wm_mask,
                                                    affine,
                                                    mtr=mtr,
                                                    ihmtr=ihmtr,
                                                    mtsat=mtsat,
                                                    ihmtsat=ihmtsat,
                                                    nufo=nufo,
                                                    bin_width=10,
                                                    fa_thr=args.fa_thr,
                                                    min_nb_voxels=args.min_nb_voxels)

    print("Saving whole brain results as txt files.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_single_fiber_mtr_ihmtr_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(w_brain_results[0], w_brain_results[1], w_brain_results[2],
                 w_brain_results[5], str(output_path), input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "WB_single_fiber_mtsat_ihmtsat_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(w_brain_results[0], w_brain_results[3], w_brain_results[4],
                 w_brain_results[5], str(output_path), input_dtype="sats")
        
    print("Fitting the whole brain results.")
    if args.in_mtr and args.in_ihmtr:
        mtr_poly_wb = fit_single_fiber_results(w_brain_results[0],
                                               w_brain_results[1],
                                               poly_order=args.poly_order)
        ihmtr_poly_wb = fit_single_fiber_results(w_brain_results[0],
                                                 w_brain_results[2],
                                                 poly_order=args.poly_order)
    if args.in_mtsat and args.in_ihmtsat:
        mtsat_poly_wb = fit_single_fiber_results(w_brain_results[0],
                                                 w_brain_results[3],
                                                 poly_order=args.poly_order)
        ihmtsat_poly_wb = fit_single_fiber_results(w_brain_results[0],
                                                   w_brain_results[4],
                                                   poly_order=args.poly_order)

    print("Saving whole brain results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_single_fiber_mtr_ihmtr_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(w_brain_results[0], w_brain_results[1], w_brain_results[2],
                   w_brain_results[5], str(output_path), mt_poly=mtr_poly_wb,
                   ihmt_poly=ihmtr_poly_wb, input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "WB_single_fiber_mtsat_ihmtsat_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(w_brain_results[0], w_brain_results[3], w_brain_results[4],
                   w_brain_results[5], str(output_path), mt_poly=mtsat_poly_wb,
                   ihmt_poly=ihmtsat_poly_wb, input_dtype="sats")

    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_single_fiber_mtr_ihmtr_bin_width_plot.png"
        output_path = out_folder / "plots" / output_name
        plot_different_bins_means([w_brain_results[0], w_brain_results_3[0], w_brain_results_5[0], w_brain_results_10[0]],
                                  [w_brain_results[1], w_brain_results_3[1], w_brain_results_5[1], w_brain_results_10[1]],
                                  [w_brain_results[2], w_brain_results_3[2], w_brain_results_5[2], w_brain_results_10[2]],
                                  [w_brain_results[5], w_brain_results_3[5], w_brain_results_5[5], w_brain_results_10[5]],
                                  str(output_path),
                                  [r'$1^\circ$ bins', r'$3^\circ$ bins', r'$5^\circ$ bins', r'$10^\circ$ bins'],
                                   input_dtype="ratios")

    print("Saving single fiber masks.")
    save_masks_by_angle_bins(e1, fa, wm_mask, affine,
                             out_folder, nufo=nufo,
                             bin_width=10, fa_thr=args.fa_thr)

#-----------------------------------CC--------------------------------------
    print("Computing single fiber averages for ROI.")
    roi_results = compute_single_fiber_averages(e1, fa,
                                                wm_mask,
                                                affine,
                                                mtr=mtr,
                                                ihmtr=ihmtr,
                                                mtsat=mtsat,
                                                ihmtsat=ihmtsat,
                                                nufo=nufo,
                                                mask=roi,
                                                bin_width=args.bin_width_mask,
                                                fa_thr=args.fa_thr,
                                                min_nb_voxels=args.min_nb_voxels)
    
    print("Saving ROI results as txt files.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "CC_single_fiber_mtr_ihmtr_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(roi_results[0], roi_results[1], roi_results[2],
                 roi_results[5], str(output_path), input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "CC_single_fiber_mtsat_ihmtsat_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(roi_results[0], roi_results[3], roi_results[4],
                 roi_results[5], str(output_path), input_dtype="sats")
        
    print("Fitting the ROI results.")
    if args.in_mtr and args.in_ihmtr:
        mtr_poly_roi = fit_single_fiber_results(roi_results[0],
                                                roi_results[1],
                                                poly_order=args.poly_order)
        ihmtr_poly_roi = fit_single_fiber_results(roi_results[0],
                                                  roi_results[2],
                                                  poly_order=args.poly_order)
    if args.in_mtsat and args.in_ihmtsat:
        mtsat_poly_roi = fit_single_fiber_results(roi_results[0],
                                                  roi_results[3],
                                                  poly_order=args.poly_order)
        ihmtsat_poly_roi = fit_single_fiber_results(roi_results[0],
                                                    roi_results[4],
                                                    poly_order=args.poly_order)

    print("Saving ROI results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "CC_single_fiber_mtr_ihmtr_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(roi_results[0], roi_results[1], roi_results[2],
                   roi_results[5], str(output_path), mt_poly=mtr_poly_roi,
                   ihmt_poly=ihmtr_poly_roi, input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "CC_single_fiber_mtsat_ihmtsat_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(roi_results[0], roi_results[3], roi_results[4],
                   roi_results[5], str(output_path), mt_poly=mtsat_poly_roi,
                   ihmt_poly=ihmtsat_poly_roi, input_dtype="sats")
        
    if args.in_mtr and args.in_ihmtr:
        output_name = "CC_&_WB_single_fiber_mtr_ihmtr_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_different_bins_means([w_brain_results[0], roi_results[0]], [w_brain_results[1], roi_results[1]], [w_brain_results[2], roi_results[2]],
                   [w_brain_results[5], roi_results[5]], str(output_path), ["WB", "CC"], mt_poly=[mtr_poly_wb, mtr_poly_roi],
                   ihmt_poly=[ihmtr_poly_wb, ihmtr_poly_roi], input_dtype="ratios")
    # if args.in_mtsat and args.in_ihmtsat:
    #     output_name = "CC_&_WB_single_fiber_mtsat_ihmtsat_plot" + files_basename + ".png"
    #     output_path = out_folder / "plots" / output_name
    #     plot_multiple_means(roi_results[0], roi_results[3], roi_results[4],
    #                roi_results[5], str(output_path), mt_poly=mtsat_poly_roi,
    #                ihmt_poly=ihmtsat_poly_roi, input_dtype="sats")

#----------------------------Crossing fibers------------------------------------
    print("Computing crossing fibers average.")
    crossing_results = compute_crossing_fibers_averages(peaks, peak_values,
                                                        wm_mask, affine,
                                                        nufo, mtr=mtr,
                                                        ihmtr=ihmtr,
                                                        mtsat=mtsat,
                                                        ihmtsat=ihmtsat,
                                                        bin_width=10,
                                                        min_nb_voxels=args.min_nb_voxels)
    
    mtr_2f_means_diag = np.diagonal(crossing_results[1], axis1=1, axis2=2)
    ihmtr_2f_means_diag = np.diagonal(crossing_results[2], axis1=1, axis2=2)
    mtsat_2f_means_diag = np.diagonal(crossing_results[3], axis1=1, axis2=2)
    ihmtsat_2f_means_diag = np.diagonal(crossing_results[4], axis1=1, axis2=2)
    nb_voxels_2f_diag = np.diagonal(crossing_results[5], axis1=1, axis2=2)

    print("Saving results as 3D plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_double_fibers_mtr_3D_plot" + files_basename
        output_path = out_folder / "plots" / output_name
        plot_3d_means(crossing_results[0], crossing_results[1][0], str(output_path), input_dtype="MTR")
        output_name = "WB_double_fibers_ihmtr_3D_plot" + files_basename
        output_path = out_folder / "plots" / output_name
        plot_3d_means(crossing_results[0], crossing_results[2][0], str(output_path), input_dtype="ihMTR")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "WB_double_fibers_mtsat_3D_plot" + files_basename
        output_path = out_folder / "plots" / output_name
        plot_3d_means(crossing_results[0], crossing_results[3][0], str(output_path), input_dtype="MTsat")
        output_name = "WB_double_fibers_ihmtsat_3D_plot" + files_basename
        output_path = out_folder / "plots" / output_name
        plot_3d_means(crossing_results[0], crossing_results[4][0], str(output_path), input_dtype="ihMTsat")
    
    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_double_fibers_mtr_ihmtr_diagonal_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_multiple_means(crossing_results[0], mtr_2f_means_diag, ihmtr_2f_means_diag,
                            nb_voxels_2f_diag, str(output_path), crossing_results[6],
                            input_dtype="ratios", legend_title=r"Peak$_1$ fraction")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "WB_double_fibers_mtsat_ihmtsat_diagonal_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_multiple_means(crossing_results[0], mtsat_2f_means_diag, ihmtsat_2f_means_diag,
                            nb_voxels_2f_diag, str(output_path), crossing_results[6],
                            input_dtype="sats", legend_title=r"Peak$_1$ fraction")
        
#----------------------------Crossing fibers analysis---------------------------
    print("Analysing crossing fibers delta_m_max.")
    if args.bin_width == 1:
        bin_min = np.array([40, 50])
        bin_max = np.array([80, 90])
        # bin_min_mid = np.mean(bin_min)
        # bin_max_mid = np.mean(bin_max)
        # min_idx = np.argwhere(w_brain_results[0] == bin_min_mid)[0][0]
        # max_idx = np.argwhere(w_brain_results[0] == bin_max_mid)[0][0]
        mtr_single_fiber_delta_m_max = np.nanmax(w_brain_results[1]) - np.nanmin(w_brain_results[1])
        # mtr_single_fiber_delta_m_max = mtr_poly_wb(bin_max_mid) - mtr_poly_wb(bin_min_mid)
        # ihmtr_single_fiber_delta_m_max = ihmtr_poly_wb(bin_min_mid) - ihmtr_poly_wb(bin_max_mid)
        # ihmtr_single_fiber_delta_m_max = w_brain_results[2][min_idx] - w_brain_results[2][max_idx]
        ihmtr_single_fiber_delta_m_max = 1

        delta_m_max_results = analyse_delta_m_max(crossing_results[0],
                                                  mtr_2f_means_diag,
                                                  ihmtr_2f_means_diag,
                                                  bin_min,
                                                  bin_max,
                                                  mtr_single_fiber_delta_m_max,
                                                  ihmtr_single_fiber_delta_m_max)
        
        mtr_poly_crossing = delta_m_max_results[0]
        # ihmtr_poly_crossing = delta_m_max_results[1]
        # mtr_poly_crossing = None
        ihmtr_poly_crossing = None
        mtsat_poly_crossing = None
        ihmtsat_poly_crossing = None
        print("Saving crossing fibers delta_m_max plot.")
        if args.in_mtr and args.in_ihmtr:
            output_name = "WB_delta_m_max_mtr_plot" + files_basename + ".png"
            output_path = out_folder / "plots" / output_name
            plot_delta_m_max(delta_m_max_results[4][1:], delta_m_max_results[2][1:],
                             delta_m_max_results[0], str(output_path))
    else:
        raise ValueError("Not implemented yet for bin_width!=1.")
    
#----------------------------3 Crossing fibers analysis-------------------------
    print("Computing 3 crossing fibers average.")
    three_crossing_results = analyse_3_crossing_fibers_averages(peaks, peak_values,
                                                        wm_mask, affine,
                                                        nufo, mtr=mtr,
                                                        ihmtr=ihmtr,
                                                        mtsat=mtsat,
                                                        ihmtsat=ihmtsat,
                                                        bin_width=30,
                                                        min_nb_voxels=args.min_nb_voxels)
    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_triple_fibers_mtr_ihmtr_diagonal_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_multiple_means(three_crossing_results[0], three_crossing_results[1], three_crossing_results[2],
                            three_crossing_results[5], str(output_path), three_crossing_results[6],
                            input_dtype="ratios", legend_title=r"Peak$_1$ fraction")

#------------------------------Correction whole brain---------------------------
    print("Correcting whole brain measures.")
    if args.in_mtr:
        corrected_mtr = correct_measure(peaks, peak_values, mtr, affine,
                                        wm_mask, mtr_poly_wb,
                                        peak_frac_thr=args.min_frac_thr,
                                        polynome_factor=mtr_poly_crossing)
        corrected_name = "mtr_corrected_WB.nii.gz"
        corrected_path = out_folder / "measures" / corrected_name
        nib.save(nib.Nifti1Image(corrected_mtr, affine), corrected_path)
    if args.in_ihmtr:
        corrected_ihmtr = correct_measure(peaks, peak_values, ihmtr, affine,
                                          wm_mask, ihmtr_poly_wb,
                                          peak_frac_thr=args.min_frac_thr,
                                          polynome_factor=ihmtr_poly_crossing)
        corrected_name = "ihmtr_corrected_WB.nii.gz"
        corrected_path = out_folder / "measures" / corrected_name
        nib.save(nib.Nifti1Image(corrected_ihmtr, affine), corrected_path)
    if args.in_mtsat:
        corrected_mtsat = correct_measure(peaks, peak_values, mtsat, affine,
                                          wm_mask, mtsat_poly_wb,
                                          peak_frac_thr=args.min_frac_thr,
                                          polynome_factor=mtsat_poly_crossing)
        corrected_name = "mtsat_corrected_WB.nii.gz"
        corrected_path = out_folder / "measures" / corrected_name
        nib.save(nib.Nifti1Image(corrected_mtsat, affine), corrected_path)
    if args.in_ihmtsat:
        corrected_ihmtsat = correct_measure(peaks, peak_values, ihmtsat,
                                            affine, wm_mask,
                                            ihmtsat_poly_wb,
                                            peak_frac_thr=args.min_frac_thr,
                                            polynome_factor=ihmtsat_poly_crossing)
        corrected_name = "ihmtsat_corrected_WB.nii.gz"
        corrected_path = out_folder / "measures" / corrected_name
        nib.save(nib.Nifti1Image(corrected_ihmtsat, affine), corrected_path)

    print("Computing difference maps.")
    if args.in_mtr:
        difference_mtr = corrected_mtr - mtr
        difference_name = "mtr_difference_WB.nii.gz"
        difference_path = out_folder / "measures" / difference_name
        nib.save(nib.Nifti1Image(difference_mtr, affine), difference_path)
    if args.in_ihmtr:
        difference_ihmtr = corrected_ihmtr - ihmtr
        difference_name = "ihmtr_difference_WB.nii.gz"
        difference_path = out_folder / "measures" / difference_name
        nib.save(nib.Nifti1Image(difference_ihmtr, affine), difference_path)
    if args.in_mtsat:
        difference_mtsat = corrected_mtsat - mtsat
        difference_name = "mtsat_difference_WB.nii.gz"
        difference_path = out_folder / "measures" / difference_name
        nib.save(nib.Nifti1Image(difference_mtsat, affine), difference_path)
    if args.in_ihmtsat:
        difference_ihmtsat = corrected_ihmtsat - ihmtsat
        difference_name = "ihmtsat_difference_WB.nii.gz"
        difference_path = out_folder / "measures" / difference_name
        nib.save(nib.Nifti1Image(difference_ihmtsat, affine), difference_path)

#----------------------------Whole brain----------------------------------------
    print("Computing single fiber averages on whole brain corrected data.")
    w_brain_cr_results = compute_single_fiber_averages(e1, fa,
                                                       wm_mask,
                                                       affine,
                                                       mtr=corrected_mtr,
                                                       ihmtr=corrected_ihmtr,
                                                       mtsat=corrected_mtsat,
                                                       ihmtsat=corrected_ihmtsat,
                                                       nufo=nufo,
                                                       bin_width=args.bin_width,
                                                       fa_thr=args.fa_thr,
                                                       min_nb_voxels=args.min_nb_voxels)
    
    print("Saving whole brain corrected results as txt files.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_single_fiber_corrected_mtr_ihmtr_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(w_brain_cr_results[0], w_brain_cr_results[1], w_brain_cr_results[2],
                 w_brain_cr_results[5], str(output_path), input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "WB_single_fiber_corrected_mtsat_ihmtsat_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(w_brain_cr_results[0], w_brain_cr_results[3], w_brain_cr_results[4],
                 w_brain_cr_results[5], str(output_path), input_dtype="sats")

    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_single_fiber_corrected_mtr_ihmtr_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(w_brain_results[0], w_brain_results[1], w_brain_results[2],
                   w_brain_results[5], str(output_path),
                   mt_cr_means=w_brain_cr_results[1],
                   ihmt_cr_means=w_brain_cr_results[2], mt_poly=mtr_poly_wb,
                   ihmt_poly=ihmtr_poly_wb, input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "WB_single_fiber_corrected_mtsat_ihmtsat_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(w_brain_results[0], w_brain_results[3], w_brain_results[4],
                   w_brain_results[5],
                   str(output_path), mt_cr_means=w_brain_cr_results[3],
                   ihmt_cr_means=w_brain_cr_results[4], mt_poly=mtsat_poly_wb,
                   ihmt_poly=ihmtsat_poly_wb, input_dtype="sats")
        
#----------------------------Bundles mask---------------------------------------
    for i in range(nb_bundles):
        ext = "".join(bundles_list[0].suffixes)
        bundle_name = str(bundles_list[i].name).replace(ext, "")
        print("Computing single fiber averages for " + bundle_name + " on corrected data.")
        bundle = nib.load(str(bundles_list[i])).get_fdata()
        bundle_results = compute_single_fiber_averages(e1, fa,
                                                       wm_mask,
                                                       affine,
                                                       mtr=mtr,
                                                       ihmtr=ihmtr,
                                                       mtsat=mtsat,
                                                       ihmtsat=ihmtsat,
                                                       nufo=nufo,
                                                       mask=bundle,
                                                       bin_width=args.bin_width_bundles,
                                                       fa_thr=args.fa_thr,
                                                       min_nb_voxels=args.min_nb_voxels)
        bundle_cr_results = compute_single_fiber_averages(e1, fa,
                                                          wm_mask,
                                                          affine,
                                                          mtr=corrected_mtr,
                                                          ihmtr=corrected_ihmtr,
                                                          mtsat=corrected_mtsat,
                                                          ihmtsat=corrected_ihmtsat,
                                                          nufo=nufo,
                                                          mask=bundle,
                                                          bin_width=args.bin_width_bundles,
                                                          fa_thr=args.fa_thr,
                                                          min_nb_voxels=args.min_nb_voxels)

        print("Saving " + bundle_name + " results as txt files.")
        if args.in_mtr and args.in_ihmtr:
            output_name = bundle_name + "_WB_single_fiber_mtr_ihmtr_results" + files_basename + ".txt"
            output_path = out_folder / "results_txt" / "bundles" / output_name
            save_txt(bundle_results[0], bundle_results[1], bundle_results[2],
                    bundle_results[5], str(output_path), input_dtype="ratios")
        if args.in_mtsat and args.in_ihmtsat:
            output_name = bundle_name + "_WB_single_fiber_mtsat_ihmtsat_results" + files_basename + ".txt"
            output_path = out_folder / "results_txt" / "bundles" / output_name
            save_txt(bundle_results[0], bundle_results[3], bundle_results[4],
                    bundle_results[5], str(output_path), input_dtype="sats")
            
        print("Saving " + bundle_name + " corrected results as txt files.")
        if args.in_mtr and args.in_ihmtr:
            output_name = bundle_name + "_WB_single_fiber_corrected_mtr_ihmtr_results" + files_basename + ".txt"
            output_path = out_folder / "results_txt" / "bundles" / output_name
            save_txt(bundle_cr_results[0], bundle_cr_results[1], bundle_cr_results[2],
                    bundle_cr_results[5], str(output_path), input_dtype="ratios")
        if args.in_mtsat and args.in_ihmtsat:
            output_name = bundle_name + "_WB_single_fiber_corrected_mtsat_ihmtsat_results" + files_basename + ".txt"
            output_path = out_folder / "results_txt" / "bundles" / output_name
            save_txt(bundle_cr_results[0], bundle_cr_results[3], bundle_cr_results[4],
                    bundle_cr_results[5], str(output_path), input_dtype="sats")
            
        print("Saving " + bundle_name + " corrected results as plots.")
        if args.in_mtr and args.in_ihmtr:
            output_name = bundle_name + "_WB_single_fiber_corrected_mtr_ihmtr_plot" + files_basename + ".png"
            output_path = out_folder / "plots" / "bundles" / output_name
            plot_means(bundle_results[0], bundle_results[1], bundle_results[2],
                       bundle_results[5], str(output_path),
                       mt_cr_means=bundle_cr_results[1],
                       ihmt_cr_means=bundle_cr_results[2], input_dtype="ratios")
        if args.in_mtsat and args.in_ihmtsat:
            output_name = bundle_name + "_WB_single_fiber_corrected_mtsat_ihmtsat_plot" + files_basename + ".png"
            output_path = out_folder / "plots" / "bundles" / output_name
            plot_means(bundle_results[0], bundle_results[3], bundle_results[4],
                       bundle_results[5], str(output_path), 
                       mt_cr_means=bundle_cr_results[3],
                       ihmt_cr_means=bundle_cr_results[4],input_dtype="sats")

#----------------------------Crossing fibers------------------------------------
    print("Computing crossing fibers average.")
    crossing_cr_results = compute_crossing_fibers_averages(peaks, peak_values,
                                                        wm_mask, affine,
                                                        nufo, mtr=corrected_mtr,
                                                        ihmtr=corrected_ihmtr,
                                                        mtsat=corrected_mtsat,
                                                        ihmtsat=corrected_ihmtsat,
                                                        bin_width=10,
                                                        min_nb_voxels=args.min_nb_voxels)
    
    mtr_cr_2f_means_diag = np.diagonal(crossing_cr_results[1], axis1=1, axis2=2)
    ihmtr_cr_2f_means_diag = np.diagonal(crossing_cr_results[2], axis1=1, axis2=2)
    mtsat_cr_2f_means_diag = np.diagonal(crossing_cr_results[3], axis1=1, axis2=2)
    ihmtsat_cr_2f_means_diag = np.diagonal(crossing_cr_results[4], axis1=1, axis2=2)
    nb_voxels_cr_2f_diag = np.diagonal(crossing_cr_results[5], axis1=1, axis2=2)

    print("Saving results as 3D plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_double_fibers_corrected_mtr_3D_plot" + files_basename
        output_path = out_folder / "plots" / output_name
        plot_3d_means(crossing_cr_results[0], crossing_cr_results[1][0], str(output_path), input_dtype="MTR")
        output_name = "WB_double_fibers_corrected_ihmtr_3D_plot" + files_basename
        output_path = out_folder / "plots" / output_name
        plot_3d_means(crossing_cr_results[0], crossing_cr_results[2][0], str(output_path), input_dtype="ihMTR")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "WB_double_fibers_corrected_mtsat_3D_plot" + files_basename
        output_path = out_folder / "plots" / output_name
        plot_3d_means(crossing_cr_results[0], crossing_cr_results[3][0], str(output_path), input_dtype="MTsat")
        output_name = "WB_double_fibers_corrected_ihmtsat_3D_plot" + files_basename
        output_path = out_folder / "plots" / output_name
        plot_3d_means(crossing_cr_results[0], crossing_cr_results[4][0], str(output_path), input_dtype="ihMTsat")
    
    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_double_fibers_corrected_mtr_ihmtr_diagonal_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_multiple_means(crossing_cr_results[0], mtr_2f_means_diag, ihmtr_2f_means_diag,
                   nb_voxels_cr_2f_diag, str(output_path), crossing_cr_results[6],
                   mt_cr_means=mtr_cr_2f_means_diag,
                   ihmt_cr_means=ihmtr_cr_2f_means_diag, input_dtype="ratios",
                   legend_title=r"Peak$_1$ fraction")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "WB_double_fibers_corrected_mtsat_ihmtsat_diagonal_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_multiple_means(crossing_cr_results[0], mtsat_2f_means_diag, ihmtsat_2f_means_diag,
                   nb_voxels_cr_2f_diag, str(output_path), crossing_cr_results[6],
                   mt_cr_means=mtsat_cr_2f_means_diag,
                   ihmt_cr_means=ihmtsat_cr_2f_means_diag, input_dtype="sats",
                   legend_title=r"Peak$_1$ fraction")
        
#----------------------------3 Crossing fibers analysis-------------------------
    print("Computing 3 crossing fibers average.")
    three_crossing_cr_results = analyse_3_crossing_fibers_averages(peaks, peak_values,
                                                        wm_mask, affine,
                                                        nufo, mtr=corrected_mtr,
                                                        ihmtr=corrected_ihmtr,
                                                        mtsat=corrected_mtsat,
                                                        ihmtsat=corrected_ihmtsat,
                                                        bin_width=30,
                                                        min_nb_voxels=args.min_nb_voxels)
    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "WB_triple_fibers_corrected_mtr_ihmtr_diagonal_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_multiple_means(three_crossing_results[0], three_crossing_results[1], three_crossing_results[2],
                            three_crossing_results[5], str(output_path), three_crossing_results[6],
                            input_dtype="ratios", legend_title=r"Peak$_1$ fraction",
                            mt_cr_means=three_crossing_cr_results[1], ihmt_cr_means=three_crossing_cr_results[2])

#------------------------------Correction CC---------------------------
    print("Correcting ROI measures.")
    if args.in_mtr:
        corrected_mtr = correct_measure(peaks, peak_values, mtr, affine,
                                        wm_mask, mtr_poly_roi,
                                        peak_frac_thr=args.min_frac_thr)
        corrected_name = "mtr_corrected_CC.nii.gz"
        corrected_path = out_folder / "measures" / corrected_name
        nib.save(nib.Nifti1Image(corrected_mtr, affine), corrected_path)
    if args.in_ihmtr:
        corrected_ihmtr = correct_measure(peaks, peak_values, ihmtr, affine,
                                          wm_mask, ihmtr_poly_roi,
                                          peak_frac_thr=args.min_frac_thr)
        corrected_name = "ihmtr_corrected_CC.nii.gz"
        corrected_path = out_folder / "measures" / corrected_name
        nib.save(nib.Nifti1Image(corrected_ihmtr, affine), corrected_path)
    if args.in_mtsat:
        corrected_mtsat = correct_measure(peaks, peak_values, mtsat, affine,
                                          wm_mask, mtsat_poly_roi,
                                          peak_frac_thr=args.min_frac_thr)
        corrected_name = "mtsat_corrected_CC.nii.gz"
        corrected_path = out_folder / "measures" / corrected_name
        nib.save(nib.Nifti1Image(corrected_mtsat, affine), corrected_path)
    if args.in_ihmtsat:
        corrected_ihmtsat = correct_measure(peaks, peak_values, ihmtsat,
                                            affine, wm_mask,
                                            ihmtsat_poly_roi,
                                            peak_frac_thr=args.min_frac_thr)
        corrected_name = "ihmtsat_corrected_CC.nii.gz"
        corrected_path = out_folder / "measures" / corrected_name
        nib.save(nib.Nifti1Image(corrected_ihmtsat, affine), corrected_path)

#----------------------------Whole brain----------------------------------------
    print("Computing single fiber averages on whole brain corrected data.")
    w_brain_cr_results = compute_single_fiber_averages(e1, fa,
                                                       wm_mask,
                                                       affine,
                                                       mtr=corrected_mtr,
                                                       ihmtr=corrected_ihmtr,
                                                       mtsat=corrected_mtsat,
                                                       ihmtsat=corrected_ihmtsat,
                                                       nufo=nufo,
                                                       bin_width=args.bin_width,
                                                       fa_thr=args.fa_thr,
                                                       min_nb_voxels=args.min_nb_voxels)
    
    print("Saving whole brain corrected results as txt files.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "CC_WB_single_fiber_corrected_mtr_ihmtr_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(w_brain_cr_results[0], w_brain_cr_results[1], w_brain_cr_results[2],
                 w_brain_cr_results[5], str(output_path), input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "CC_WB_single_fiber_corrected_mtsat_ihmtsat_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(w_brain_cr_results[0], w_brain_cr_results[3], w_brain_cr_results[4],
                 w_brain_cr_results[5], str(output_path), input_dtype="sats")

    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "CC_WB_single_fiber_corrected_mtr_ihmtr_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(w_brain_results[0], w_brain_results[1], w_brain_results[2],
                   w_brain_results[5], str(output_path),
                   mt_cr_means=w_brain_cr_results[1],
                   ihmt_cr_means=w_brain_cr_results[2], mt_poly=mtr_poly_wb,
                   ihmt_poly=ihmtr_poly_wb, input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "CC_WB_single_fiber_corrected_mtsat_ihmtsat_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(w_brain_results[0], w_brain_results[3], w_brain_results[4],
                   w_brain_results[5],
                   str(output_path), mt_cr_means=w_brain_cr_results[3],
                   ihmt_cr_means=w_brain_cr_results[4], mt_poly=mtsat_poly_wb,
                   ihmt_poly=ihmtsat_poly_wb, input_dtype="sats")
        
#-----------------------------------CC----------------------------------------
    print("Computing single fiber averages on ROI corrected data.")
    roi_cr_results = compute_single_fiber_averages(e1, fa,
                                                   wm_mask,
                                                   affine,
                                                   mtr=corrected_mtr,
                                                   ihmtr=corrected_ihmtr,
                                                   mtsat=corrected_mtsat,
                                                   ihmtsat=corrected_ihmtsat,
                                                   nufo=nufo,
                                                   mask=roi,
                                                   bin_width=args.bin_width_mask,
                                                   fa_thr=args.fa_thr,
                                                   min_nb_voxels=args.min_nb_voxels)
    
    print("Saving ROI corrected results as txt files.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "CC_single_fiber_corrected_mtr_ihmtr_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(roi_cr_results[0], roi_cr_results[1], roi_cr_results[2],
                 roi_cr_results[5], str(output_path), input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "CC_single_fiber_corrected_mtsat_ihmtsat_results" + files_basename + ".txt"
        output_path = out_folder / "results_txt" / output_name
        save_txt(roi_cr_results[0], roi_cr_results[3], roi_cr_results[4],
                 roi_cr_results[5], str(output_path), input_dtype="sats")

    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "CC_single_fiber_corrected_mtr_ihmtr_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(roi_results[0], roi_results[1], roi_results[2],
                   roi_results[5], str(output_path),
                   mt_cr_means=roi_cr_results[1],
                   ihmt_cr_means=roi_cr_results[2], mt_poly=mtr_poly_roi,
                   ihmt_poly=ihmtr_poly_roi, input_dtype="ratios")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "CC_single_fiber_corrected_mtsat_ihmtsat_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_means(roi_results[0], roi_results[3], roi_results[4],
                   roi_results[5],
                   str(output_path), mt_cr_means=roi_cr_results[3],
                   ihmt_cr_means=roi_cr_results[4], mt_poly=mtsat_poly_roi,
                   ihmt_poly=ihmtsat_poly_roi, input_dtype="sats")
        
#----------------------------Bundles mask---------------------------------------
    for i in range(nb_bundles):
        ext = "".join(bundles_list[0].suffixes)
        bundle_name = str(bundles_list[i].name).replace(ext, "")
        print("Computing single fiber averages for " + bundle_name + " on corrected data.")
        bundle = nib.load(str(bundles_list[i])).get_fdata()
        bundle_results = compute_single_fiber_averages(e1, fa,
                                                       wm_mask,
                                                       affine,
                                                       mtr=mtr,
                                                       ihmtr=ihmtr,
                                                       mtsat=mtsat,
                                                       ihmtsat=ihmtsat,
                                                       nufo=nufo,
                                                       mask=bundle,
                                                       bin_width=args.bin_width_bundles,
                                                       fa_thr=args.fa_thr,
                                                       min_nb_voxels=args.min_nb_voxels)
        bundle_cr_results = compute_single_fiber_averages(e1, fa,
                                                          wm_mask,
                                                          affine,
                                                          mtr=corrected_mtr,
                                                          ihmtr=corrected_ihmtr,
                                                          mtsat=corrected_mtsat,
                                                          ihmtsat=corrected_ihmtsat,
                                                          nufo=nufo,
                                                          mask=bundle,
                                                          bin_width=args.bin_width_bundles,
                                                          fa_thr=args.fa_thr,
                                                          min_nb_voxels=args.min_nb_voxels)

        print("Saving " + bundle_name + " results as txt files.")
        if args.in_mtr and args.in_ihmtr:
            output_name = bundle_name + "_CC_single_fiber_mtr_ihmtr_results" + files_basename + ".txt"
            output_path = out_folder / "results_txt" / "bundles" / output_name
            save_txt(bundle_results[0], bundle_results[1], bundle_results[2],
                    bundle_results[5], str(output_path), input_dtype="ratios")
        if args.in_mtsat and args.in_ihmtsat:
            output_name = bundle_name + "_CC_single_fiber_mtsat_ihmtsat_results" + files_basename + ".txt"
            output_path = out_folder / "results_txt" / "bundles" / output_name
            save_txt(bundle_results[0], bundle_results[3], bundle_results[4],
                    bundle_results[5], str(output_path), input_dtype="sats")
            
        print("Saving " + bundle_name + " corrected results as txt files.")
        if args.in_mtr and args.in_ihmtr:
            output_name = bundle_name + "_CC_single_fiber_corrected_mtr_ihmtr_results" + files_basename + ".txt"
            output_path = out_folder / "results_txt" / "bundles" / output_name
            save_txt(bundle_cr_results[0], bundle_cr_results[1], bundle_cr_results[2],
                    bundle_cr_results[5], str(output_path), input_dtype="ratios")
        if args.in_mtsat and args.in_ihmtsat:
            output_name = bundle_name + "_CC_single_fiber_corrected_mtsat_ihmtsat_results" + files_basename + ".txt"
            output_path = out_folder / "results_txt" / "bundles" / output_name
            save_txt(bundle_cr_results[0], bundle_cr_results[3], bundle_cr_results[4],
                    bundle_cr_results[5], str(output_path), input_dtype="sats")
            
        print("Saving " + bundle_name + " corrected results as plots.")
        if args.in_mtr and args.in_ihmtr:
            output_name = bundle_name + "_CC_single_fiber_corrected_mtr_ihmtr_plot" + files_basename + ".png"
            output_path = out_folder / "plots" / "bundles" / output_name
            plot_means(bundle_results[0], bundle_results[1], bundle_results[2],
                       bundle_results[5], str(output_path),
                       mt_cr_means=bundle_cr_results[1],
                       ihmt_cr_means=bundle_cr_results[2], input_dtype="ratios")
        if args.in_mtsat and args.in_ihmtsat:
            output_name = bundle_name + "_CC_single_fiber_corrected_mtsat_ihmtsat_plot" + files_basename + ".png"
            output_path = out_folder / "plots" / "bundles" / output_name
            plot_means(bundle_results[0], bundle_results[3], bundle_results[4],
                       bundle_results[5], str(output_path), 
                       mt_cr_means=bundle_cr_results[3],
                       ihmt_cr_means=bundle_cr_results[4],input_dtype="sats")

#----------------------------Crossing fibers------------------------------------
    print("Computing crossing fibers average.")
    crossing_cr_results = compute_crossing_fibers_averages(peaks, peak_values,
                                                        wm_mask, affine,
                                                        nufo, mtr=corrected_mtr,
                                                        ihmtr=corrected_ihmtr,
                                                        mtsat=corrected_mtsat,
                                                        ihmtsat=corrected_ihmtsat,
                                                        bin_width=10,
                                                        min_nb_voxels=args.min_nb_voxels)
    
    mtr_cr_2f_means_diag = np.diagonal(crossing_cr_results[1], axis1=1, axis2=2)
    ihmtr_cr_2f_means_diag = np.diagonal(crossing_cr_results[2], axis1=1, axis2=2)
    mtsat_cr_2f_means_diag = np.diagonal(crossing_cr_results[3], axis1=1, axis2=2)
    ihmtsat_cr_2f_means_diag = np.diagonal(crossing_cr_results[4], axis1=1, axis2=2)
    nb_voxels_cr_2f_diag = np.diagonal(crossing_cr_results[5], axis1=1, axis2=2)
    
    print("Saving results as plots.")
    if args.in_mtr and args.in_ihmtr:
        output_name = "CC_double_fibers_corrected_mtr_ihmtr_diagonal_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_multiple_means(crossing_cr_results[0], mtr_2f_means_diag, ihmtr_2f_means_diag,
                   nb_voxels_cr_2f_diag, str(output_path), crossing_cr_results[6],
                   mt_cr_means=mtr_cr_2f_means_diag,
                   ihmt_cr_means=ihmtr_cr_2f_means_diag, input_dtype="ratios",
                   legend_title=r"Peak$_1$ fraction")
    if args.in_mtsat and args.in_ihmtsat:
        output_name = "CC_double_fibers_corrected_mtsat_ihmtsat_diagonal_plot" + files_basename + ".png"
        output_path = out_folder / "plots" / output_name
        plot_multiple_means(crossing_cr_results[0], mtsat_2f_means_diag, ihmtsat_2f_means_diag,
                   nb_voxels_cr_2f_diag, str(output_path), crossing_cr_results[6],
                   mt_cr_means=mtsat_cr_2f_means_diag,
                   ihmt_cr_means=ihmtsat_cr_2f_means_diag, input_dtype="sats",
                   legend_title=r"Peak$_1$ fraction")

    # print("Computing poor's man ihMTR.")
    # if args.in_mtr:
    #     ihmtr_poor = compute_poor_ihmtr(corrected_mtr, wm_mask, mtr_poly)
    #     ihmtr_poor_name = "poor_ihmtr.nii.gz"
    #     ihmtr_poor_path = out_folder / ihmtr_poor_name
    #     nib.save(nib.Nifti1Image(ihmtr_poor, affine), ihmtr_poor_path)

    # print("Saving masked versions.")
    # wm_mask_bool = (wm_mask > 0.9)
    # if args.in_mtr:
    #     masked_corrected_mtr = np.where(wm_mask_bool, corrected_mtr, 0)
    #     corrected_name = "mtr_corrected_masked.nii.gz"
    #     corrected_path = out_folder / corrected_name
    #     nib.save(nib.Nifti1Image(masked_corrected_mtr, affine), corrected_path)
    #     masked_mtr = np.where(wm_mask_bool, mtr, 0)
    #     corrected_name = "mtr_masked.nii.gz"
    #     corrected_path = out_folder / corrected_name
    #     nib.save(nib.Nifti1Image(masked_mtr, affine), corrected_path)
    # if args.in_ihmtr:
    #     masked_corrected_ihmtr = np.where(wm_mask_bool, corrected_ihmtr, 0)
    #     corrected_name = "ihmtr_corrected_masked.nii.gz"
    #     corrected_path = out_folder / corrected_name
    #     nib.save(nib.Nifti1Image(masked_corrected_ihmtr, affine), corrected_path)
    #     masked_ihmtr = np.where(wm_mask_bool, ihmtr, 0)
    #     corrected_name = "ihmtr_masked.nii.gz"
    #     corrected_path = out_folder / corrected_name
    #     nib.save(nib.Nifti1Image(masked_ihmtr, affine), corrected_path)
    # if args.in_mtsat:
    #     masked_corrected_mtsat = np.where(wm_mask_bool, corrected_mtsat, 0)
    #     corrected_name = "mtsat_corrected_masked.nii.gz"
    #     corrected_path = out_folder / corrected_name
    #     nib.save(nib.Nifti1Image(masked_corrected_mtsat, affine), corrected_path)
    #     masked_mtsat = np.where(wm_mask_bool, mtsat, 0)
    #     corrected_name = "mtsat_masked.nii.gz"
    #     corrected_path = out_folder / corrected_name
    #     nib.save(nib.Nifti1Image(masked_mtsat, affine), corrected_path)
    # if args.in_ihmtsat:
    #     masked_corrected_ihmtsat = np.where(wm_mask_bool, corrected_ihmtsat, 0)
    #     corrected_name = "ihmtsat_corrected_masked.nii.gz"
    #     corrected_path = out_folder / corrected_name
    #     nib.save(nib.Nifti1Image(masked_corrected_ihmtsat, affine), corrected_path)
    #     masked_ihmtsat = np.where(wm_mask_bool, ihmtsat, 0)
    #     corrected_name = "ihmtsat_masked.nii.gz"
    #     corrected_path = out_folder / corrected_name
    #     nib.save(nib.Nifti1Image(masked_ihmtsat, affine), corrected_path)

if __name__ == "__main__":
    main()
