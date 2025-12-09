#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to plot the track-profile of a bundle for both the MTR and the
fixel-wise MTR.
"""

import argparse
import logging
import warnings

import scipy.stats as stats
from matplotlib.gridspec import GridSpec
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D 
import matplotlib as mpl
from matplotlib import pyplot as plt
import nibabel as nib
import numpy as np

from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
     p = argparse.ArgumentParser(
          description=__doc__,
          formatter_class=argparse.RawTextHelpFormatter)

     p.add_argument('in_bundle_name', type=str,
                    help='Name of the bundle being processed.')

     p.add_argument('out_dir',
                    help='Output directory to save the track-profiles plot.')

     # Inputs for all subjects (scan versions)
     p.add_argument('--in_mtr_profiles_all', nargs='+', required=True,
                    help='Input bundle-specific MTR track-profile txt files '
                         'for all subjects (scan).')

     p.add_argument('--in_fixel_mtr_profiles_all', nargs='+', required=True,
                    help='Input bundle-specific fixel-wise MTR track-profile '
                         'txt files for all subjects (scan).')

     p.add_argument('--in_afd_fixel_profiles_all', nargs='+', required=True,
                    help='Input bundle-specific AFD fixel track-profile txt '
                         'files for all subjects (scan).')

     p.add_argument('--in_nufo_profiles_all', nargs='+', required=True,
                    help='Input NuFO track-profile txt files for all subjects '
                         '(scan).')

     p.add_argument('--in_nb_voxels_profiles_all', nargs='+',
                    help='Input NuFO track-profile txt files for all subjects '
                         '(scan). OPTIONAL.')

     # Inputs for subjects with rescans (scan versions)
     p.add_argument('--in_mtr_profiles_scan', nargs='+', required=True,
                    help='Input bundle-specific MTR track-profile txt files '
                         'for the scan version of subjects with rescans.')

     p.add_argument('--in_fixel_mtr_profiles_scan', nargs='+', required=True,
                    help='Input bundle-specific fixel-wise MTR track-profile '
                         'txt files for the scan version of subjects with '
                         'rescans')

     # Inputs for subjects with rescans (rescan versions)
     p.add_argument('--in_mtr_profiles_rescan', nargs='+', required=True,
                    help='Input bundle-specific MTR track-profile txt files '
                         'for the rescan version of subjects with rescans.')

     p.add_argument('--in_fixel_mtr_profiles_rescan', nargs='+', required=True,
                    help='Input bundle-specific fixel-wise MTR track-profile '
                         'txt files for the rescan version of subjects with '
                         'rescans')

     # Other parameters
     p.add_argument('--nb_sections', type=int, default=20,
                    help='Number of sections in the bundle labels. '
                         'Default is 20.')
    
     p.add_argument('--min_nb_subjects', type=int, default=5,
                    help='Minimum number of subjects required to compute the '
                         'mean profile value at each bundle section. '
                         'Default is 5.')
     
     # Add-ons
     p.add_argument('--in_overlap_txt', type=str,
                    help='TXT file with bundle crossings > threshold.')

     p.add_argument('--in_significance_txt', type=str,
                    help='TXT file with significant sections per bundle.')
     
     p.add_argument('--in_mtr_profiles_overlap', nargs='+', action='append',
                    help='Input bundle-specific MTR track-profile txt files '
                         'from the overlap analysis.')

     p.add_argument('--in_fixel_mtr_profiles_overlap', nargs='+',
                    action='append',
                    help='Input bundle-specific fixel-wise MTR track-profile '
                         'txt files from the overlap analysis.')
     
     p.add_argument('--in_nb_voxels_profiles_overlap', nargs='+',
                    action='append',
                    help='Input NuFO track-profile txt files from the overlap '
                         'analysis.')

     add_verbose_arg(p)
     add_overwrite_arg(p)

     return p


def load_overlap_sections_from_txt(txt_file, target_bundle):
     overlap_sections = set()
     if txt_file is None:
          return overlap_sections

     with open(txt_file, "r") as f:
          for line in f:
               line = line.strip()
               # Skip headers
               if not line or line.startswith("#") or "-->" not in line:
                    continue
               left, _= line.split("-->")
               # Left side: "BUNDLE - section X"
               bundle_part = left.split("- section")[0].strip()
               section_part = left.split("- section")[1].strip()
               if bundle_part != target_bundle:
                    continue
               sec = int(section_part)
               overlap_sections.add(sec)

     return overlap_sections


def load_significant_sections_from_txt(txt_file, target_bundle):
     significant_sections = set()
     if txt_file is None:
          return significant_sections

     with open(txt_file, "r") as f:
          lines = f.readlines()
     inside_bundle_block = False

     for line in lines:
          line = line.strip()
          if line.startswith("==="):
               inside_bundle_block = target_bundle in line
               continue
          if inside_bundle_block:
               if "No significant sections" in line:
                    break
               if line.startswith("Section"):
                    sec = int(line.split("Section")[1].strip())
                    significant_sections.add(sec)

     return significant_sections


def main():
     parser = _build_arg_parser()
     args = parser.parse_args()

     if args.verbose:
          logging.basicConfig(level=logging.INFO)

     warnings.filterwarnings('ignore')

     assert_inputs_exist(parser, [args.in_mtr_profiles_all,
                                  args.in_fixel_mtr_profiles_all,
                                  args.in_afd_fixel_profiles_all,
                                  args.in_nufo_profiles_all,
                                  args.in_mtr_profiles_scan,
                                  args.in_fixel_mtr_profiles_scan,
                                  args.in_mtr_profiles_rescan,
                                  args.in_fixel_mtr_profiles_rescan])
     assert_outputs_exist(parser, args, [args.out_dir])

     # Load profiles
     mtr_profiles_all = np.array([np.loadtxt(f) for f in args.in_mtr_profiles_all])
     fixel_mtr_profiles_all = np.array([np.loadtxt(f) for f in args.in_fixel_mtr_profiles_all])
     nufo_profiles_all = np.array([np.loadtxt(f) for f in args.in_nufo_profiles_all])
     afd_profiles_all = np.array([np.loadtxt(f) for f in args.in_afd_fixel_profiles_all])
     mtr_profiles_scan = np.array([np.loadtxt(f) for f in args.in_mtr_profiles_scan])
     fixel_mtr_profiles_scan = np.array([np.loadtxt(f) for f in args.in_fixel_mtr_profiles_scan])
     mtr_profiles_rescan = np.array([np.loadtxt(f) for f in args.in_mtr_profiles_rescan])
     fixel_mtr_profiles_rescan = np.array([np.loadtxt(f) for f in args.in_fixel_mtr_profiles_rescan])
     if args.in_nb_voxels_profiles_all:
          nb_voxels_profiles_all = np.array([np.loadtxt(f) for f in args.in_nb_voxels_profiles_all])
     labels = np.arange(1, args.nb_sections + 1, 1)

     # Prepare principal profiles by averaging all scans
     nb_subjects = mtr_profiles_all.shape[0]
     mtr_profiles = np.where(mtr_profiles_all > 0, mtr_profiles_all, np.nan)
     fixel_mtr_profiles = np.where(fixel_mtr_profiles_all > 0,
                                   fixel_mtr_profiles_all, np.nan)
     nufo_profiles = np.where(nufo_profiles_all > 0, nufo_profiles_all, np.nan)
     afd_profiles = np.where(afd_profiles_all > 0, afd_profiles_all, np.nan)
     # Compute number of subjects with valid data at each section
     nb_subjects_all = np.sum(~np.isnan(fixel_mtr_profiles), axis=0)
     min_nb_subjects_mask = nb_subjects_all >= args.min_nb_subjects
     # Compute mean and std profiles
     mtr_profile = np.nanmean(mtr_profiles, axis=0)
     fixel_mtr_profile = np.nanmean(fixel_mtr_profiles, axis=0)
     nufo_profile = np.nanmean(nufo_profiles, axis=0)
     afd_profile = np.nanmean(afd_profiles, axis=0)
     mtr_profile_std = np.nanstd(mtr_profiles, axis=0)
     fixel_mtr_profile_std = np.nanstd(fixel_mtr_profiles, axis=0)
     mtr_profile = np.where(np.isnan(mtr_profile), 0, mtr_profile)
     fixel_mtr_profile = np.where(np.isnan(fixel_mtr_profile), 0,
                                  fixel_mtr_profile)
     nufo_profile = np.where(np.isnan(nufo_profile), 0, nufo_profile)
     afd_profile = np.where(np.isnan(afd_profile), 0, afd_profile)
     if args.in_nb_voxels_profiles_all:
          nb_voxels_profiles = np.where(nb_voxels_profiles_all > 0,
                                        nb_voxels_profiles_all, np.nan)
          nb_voxels_profile = np.nanmean(nb_voxels_profiles, axis=0)
          nb_voxels_profile = np.where(np.isnan(nb_voxels_profile), 0,
                                       nb_voxels_profile)

     # Compute the absolute difference between scan and rescan profiles
     nb_subjects_scan = mtr_profiles_scan.shape[0]
     mtr_profile_diff = np.zeros((nb_subjects_scan, args.nb_sections))
     fixel_mtr_profile_diff = np.zeros((nb_subjects_scan, args.nb_sections))
     for i in range(nb_subjects_scan):
          mtr_profile_diff[i, :] = np.abs(mtr_profiles_scan[i] - mtr_profiles_rescan[i]) / mtr_profiles_scan[i] * 100
          fixel_mtr_profile_diff[i, :] = np.abs(fixel_mtr_profiles_scan[i] - fixel_mtr_profiles_rescan[i]) / fixel_mtr_profiles_scan[i] * 100
     # Masks to consider only valid data points
     mtr_mask = (mtr_profiles_scan != 0) & (mtr_profiles_rescan != 0) & (mtr_profile_diff != np.inf) & (~np.isnan(mtr_profile_diff))
     fixel_mtr_mask = (fixel_mtr_profiles_scan != 0) & (fixel_mtr_profiles_rescan != 0) & (fixel_mtr_profile_diff != np.inf) & (~np.isnan(fixel_mtr_profile_diff))
     # Masks to consider only sections with minimum data points (subjects)
     mtr_mask = mtr_mask & np.repeat((np.sum(mtr_mask, axis=0) >= args.min_nb_subjects)[np.newaxis, :], nb_subjects_scan, axis=0)
     fixel_mtr_mask = fixel_mtr_mask & np.repeat((np.sum(fixel_mtr_mask, axis=0) >= args.min_nb_subjects)[np.newaxis, :], nb_subjects_scan, axis=0)

     # Load overlap and significant sections
     overlap_sections = load_overlap_sections_from_txt(
     args.in_overlap_txt, args.in_bundle_name)
     significant_sections = load_significant_sections_from_txt(
     args.in_significance_txt, args.in_bundle_name)

     # Load overlap profiles if provided
     if (args.in_mtr_profiles_overlap is not None) and (args.in_fixel_mtr_profiles_overlap is not None):
          # Load MTR
          mtr_profiles_overlap = np.zeros((len(args.in_mtr_profiles_overlap),
                                           nb_subjects, args.nb_sections))
          for i, overlap_profiles in enumerate(args.in_mtr_profiles_overlap):
               mtr_profiles_overlap[i] = np.array([np.loadtxt(f) for f in overlap_profiles])
          # Load fixel-wise MTR
          fixel_mtr_profiles_overlap = np.zeros((len(args.in_fixel_mtr_profiles_overlap),
                                                 nb_subjects, args.nb_sections))
          for i, overlap_profiles in enumerate(args.in_fixel_mtr_profiles_overlap):
               fixel_mtr_profiles_overlap[i] = np.array([np.loadtxt(f) for f in overlap_profiles])
          # Process overlap profiles
          mtr_profiles_overlap = np.where(mtr_profiles_overlap > 0,
                                          mtr_profiles_overlap, np.nan)
          fixel_mtr_profiles_overlap = np.where(fixel_mtr_profiles_overlap > 0,
                                                fixel_mtr_profiles_overlap,
                                                np.nan)
          # ttest_overlap = stats.ttest_rel(mtr_profiles_overlap,
          #                                 fixel_mtr_profiles_overlap,
          #                                 axis=1, nan_policy='omit')
          # pvalues_overlap = stats.false_discovery_control(ttest_overlap.pvalue,
          #                                                 axis=1, method='bh')
          nb_subjects_overlap = np.sum(~np.isnan(fixel_mtr_profiles_overlap),
                                       axis=1)
          min_nb_subjects_overlap_mask = nb_subjects_overlap >= args.min_nb_subjects
          mtr_profiles_overlap = np.nanmean(mtr_profiles_overlap, axis=1)
          fixel_mtr_profiles_overlap = np.nanmean(fixel_mtr_profiles_overlap,
                                                  axis=1)
          mtr_profiles_overlap = np.where(np.isnan(mtr_profiles_overlap), 0,
                                          mtr_profiles_overlap)
          fixel_mtr_profiles_overlap = np.where(np.isnan(fixel_mtr_profiles_overlap),
                                                0, fixel_mtr_profiles_overlap)
     if args.in_nb_voxels_profiles_overlap is not None:
          nb_voxels_profiles_overlap = np.zeros((len(args.in_nb_voxels_profiles_overlap),
                                                 nb_subjects, args.nb_sections))
          for i, overlap_profiles in enumerate(args.in_nb_voxels_profiles_overlap):
               nb_voxels_profiles_overlap[i] = np.array([np.loadtxt(f) for f in overlap_profiles])
          nb_voxels_profiles_overlap = np.where(nb_voxels_profiles_overlap > 0,
                                                nb_voxels_profiles_overlap,
                                                np.nan)
          nb_voxels_profiles_overlap = np.nanmean(nb_voxels_profiles_overlap, axis=1)
          nb_voxels_profiles_overlap = np.where(np.isnan(nb_voxels_profiles_overlap), 0,
                                                nb_voxels_profiles_overlap)

     data_for_boxplot = []
     positions = []
     section_centers = []
     offset = 0.15  # distance between MTR and Fixel boxes (smaller = closer)
     for sec in range(args.nb_sections):
          # MTR diff values at this section across subjects
          data_for_boxplot.append(mtr_profile_diff[mtr_mask[:, sec], sec])
          # Fixel-MTR diff
          data_for_boxplot.append(fixel_mtr_profile_diff[fixel_mtr_mask[:, sec], sec])
          center = sec + 1  # section index on x-axis
          section_centers.append(center)
          positions.append(center - offset)  # MTR
          positions.append(center + offset)  # Fixel-MTR

     # Plot profiles
     colors = ['#00A759', '#B45E2F']
     norm = mpl.colors.Normalize(vmin=0.3, vmax=0.7)

     fig = plt.figure(figsize=(8, 7))
     gs = GridSpec(2, 2, width_ratios=[20, 1], height_ratios=[3, 1],
                   wspace=0.15, hspace=0.25)
     
     ax1 = fig.add_subplot(gs[0, 0])
     ymin = 0.25
     ymax = 0.55

     ax1.fill_between(labels[min_nb_subjects_mask],
                      mtr_profile[min_nb_subjects_mask] - mtr_profile_std[min_nb_subjects_mask],
                      mtr_profile[min_nb_subjects_mask] + mtr_profile_std[min_nb_subjects_mask],
                      color=colors[0], alpha=0.2)
     ax1.fill_between(labels[min_nb_subjects_mask],
                      fixel_mtr_profile[min_nb_subjects_mask] - fixel_mtr_profile_std[min_nb_subjects_mask],
                      fixel_mtr_profile[min_nb_subjects_mask] + fixel_mtr_profile_std[min_nb_subjects_mask],
                      color=colors[1], alpha=0.2)

     # --- Plot with markers ---
     ax1.plot(labels[min_nb_subjects_mask], mtr_profile[min_nb_subjects_mask],
              label='MTR', color=colors[0])
     ax1.plot(labels[min_nb_subjects_mask],
              fixel_mtr_profile[min_nb_subjects_mask], label='Fixel-wise MTR',
              color=colors[1])
     for sec in range(1, args.nb_sections + 1):
          if not min_nb_subjects_mask[sec - 1]:
               continue
          x = sec
          y_mtr = mtr_profile[sec - 1]
          y_fixel = fixel_mtr_profile[sec - 1]
          # Use 'x' if overlap, else 'o'
          marker_style = 'X' if sec in overlap_sections else 'o'
          ax1.scatter(x, y_mtr, marker=marker_style, color=colors[0], zorder=5)
          ax1.scatter(x, y_fixel, marker=marker_style, color=colors[1],
                      zorder=5)

     # Double-sided arrows for significance
     linewidth = 1.0
     cap_half_width = 0.05   # half-width of horizontal cap in X units
     for sec in sorted(significant_sections):
          if not min_nb_subjects_mask[sec - 1]:
               continue
          x = sec
          y1 = mtr_profile[sec - 1]
          y2 = fixel_mtr_profile[sec - 1]
          # ---- Order endpoints
          y_low, y_high = sorted([y1, y2])
          # ---- Clip to visible Y range
          y_low_clip = max(y_low, ymin)
          y_high_clip = min(y_high, ymax)
          # Fully invisible -> skip
          if y_low_clip >= y_high_clip:
               continue
          # ---- Central vertical line (always clipped)
          ax1.plot([x, x], [y_low_clip, y_high_clip], color='black',
                   linewidth=linewidth, zorder=6, clip_on=True)
          # ---- LOWER FLAT CAP (only if true endpoint is inside)
          if y_low >= ymin:
               ax1.plot([x - cap_half_width, x + cap_half_width],
                        [y_low, y_low], color='black', linewidth=linewidth,
                        zorder=7, clip_on=True)
          # ---- UPPER FLAT CAP (only if true endpoint is inside)
          if y_high <= ymax:
               ax1.plot([x - cap_half_width, x + cap_half_width],
                        [y_high, y_high], color='black', linewidth=linewidth,
                        zorder=7, clip_on=True)

     # Add secondary axis for NuFO
     ax2 = ax1.twinx()
     ax2.plot(labels[min_nb_subjects_mask],
              nufo_profile[min_nb_subjects_mask],
              linestyle='--', color="darkgrey", zorder=3)
     ax2.scatter(labels[min_nb_subjects_mask],
                 nufo_profile[min_nb_subjects_mask],
                 c=afd_profile[min_nb_subjects_mask],
                 cmap='Greys', norm=norm, edgecolors='darkgrey', zorder=4)
     ax2.set_ylabel('Mean NuFO', color="darkgrey")
     ax2.tick_params(axis='y', labelcolor="darkgrey")

     # Add AFD colorbar
     ax_cb = fig.add_subplot(gs[0, 1])
     sm = plt.cm.ScalarMappable(cmap='Greys', norm=norm)
     sm.set_array([])
     cbar = fig.colorbar(sm, cax=ax_cb)
     cbar.set_label('Mean AFD', color="darkgrey")
     cbar.ax.tick_params(color='darkgrey', labelcolor='darkgrey')
     cbar.outline.set_edgecolor('darkgrey')
     # cbar.tick_params(axis='y', labelcolor="darkgrey")

     # Plot overlap markers if provided
     if (args.in_mtr_profiles_overlap is not None) and (args.in_fixel_mtr_profiles_overlap is not None):
          markers_overlap = ['x', '*', '+', 'D', 's']
          for i in range(mtr_profiles_overlap.shape[0]):
               ax1.scatter(labels[min_nb_subjects_overlap_mask[i]],
                           mtr_profiles_overlap[i][min_nb_subjects_overlap_mask[i]],
                           marker=markers_overlap[i],
                           color=colors[0], zorder=5)
               ax1.scatter(labels[min_nb_subjects_overlap_mask[i]],
                           fixel_mtr_profiles_overlap[i][min_nb_subjects_overlap_mask[i]],
                           marker=markers_overlap[i],
                           color=colors[1], zorder=5)

     # Axis labels
     ax1.set_xlabel('Bundle section')
     ax1.set_ylabel('Mean MTR')
     ax1.set_title('Track-profile of MTR and fixel-wise MTR for the {} bundle'.format(args.in_bundle_name))

     # Set legend
     legend_handles = []
     # ---- Main profiles (circle + line) ----
     legend_handles.append(Line2D([0], [0], color=colors[0], marker='o',
                                  linestyle='-', label='MTR'))
     legend_handles.append(Line2D([0], [0], color=colors[1], marker='o',
                                  linestyle='-', label='Fixel-wise MTR'))

     # ---- Complexity (NuFO + AFD color) ----
     legend_handles.append(Line2D([0], [0], color='darkgrey', linestyle='--',
                                  marker='o', markerfacecolor='white',
                                  markeredgecolor='darkgrey',
                                  label='Complexity'))

     # ---- Overlap as "X" marker ----
     if args.in_overlap_txt:
          overlap_handle = (Line2D([0], [0], color=colors[0], marker='X',
                                   linestyle='None'),
                            Line2D([0], [0], color=colors[1], marker='X',
                                   linestyle='None'))
          legend_handles.append(overlap_handle)

     ax1.legend(handles=legend_handles,
                labels=[h.get_label() if not isinstance(h, tuple) else ">20% overlap" for h in legend_handles],
                handler_map={tuple: HandlerTuple(ndivide=None)},
                loc='best')

     ax1.set_ylim(ymin, ymax)
     ax1.set_xlim(0, args.nb_sections + 1)
     ax1.set_xticks(np.arange(1, args.nb_sections + 1, 1))
     ax2.set_ylim(1, 5)
     ax2.set_yticks(np.arange(1, 6, 1))

     # Absolute difference subplot
     ax3 = fig.add_subplot(gs[1, 0])

     # Create boxplot
     bp = ax3.boxplot(data_for_boxplot, positions=positions, patch_artist=True,
                      showfliers=False, medianprops=dict(color='black'), widths=0.2)
     # Color the boxes alternating (MTR = cmap[0], Fixel = cmap[1])
     for i, box in enumerate(bp['boxes']):
          col = colors[0] if i % 2 == 0 else colors[1]
          box.set_facecolor(col)
          box.set_alpha(0.5)

     ax3.set_xticks(section_centers)
     ax3.set_xticklabels(section_centers)

     ax3.set_ylabel('|%diff| scan-rescan')
     ax3.set_xlabel('Bundle section')
     ax3.set_ylim(0, 20)
     ax3.set_yticks([0, 5, 10, 15, 20])
     ax3.set_xlim(0, args.nb_sections + 1)
     # ax3.set_xticks(np.arange(1, args.nb_sections + 1, 1))
     # ax3.legend(loc='upper right')
     # ax3.grid(alpha=0.3)

     plt.tight_layout()
     # plt.show()
     plt.savefig(args.out_dir + '/track_profile_{}.png'.format(args.in_bundle_name),
                 dpi=500)


if __name__ == "__main__":
    main()
