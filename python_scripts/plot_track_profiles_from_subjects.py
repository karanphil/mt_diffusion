#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to plot the track-profile of a bundle for both the MTR and the
fixel-wise MTR.
"""

import argparse
import logging
import warnings

from matplotlib.gridspec import GridSpec
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

     add_verbose_arg(p)
     add_overwrite_arg(p)

     return p


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
                                  args.in_fixel_mtr_profiles_rescan,])
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
     labels = np.arange(1, args.nb_sections + 1, 1)

     # Prepare principal profiles by averaging all scans
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

     # Compute the absolute difference between scan and rescan profiles
     nb_subjects = mtr_profiles_scan.shape[0]
     mtr_profile_diff = np.zeros((nb_subjects, args.nb_sections))
     fixel_mtr_profile_diff = np.zeros((nb_subjects, args.nb_sections))
     for i in range(nb_subjects):
          # Might have to add a check to avoid comparing empty sections
          mtr_profile_diff[i, :] = np.abs(mtr_profiles_scan[i] - mtr_profiles_rescan[i]) / mtr_profiles_scan[i] * 100
          fixel_mtr_profile_diff[i, :] = np.abs(fixel_mtr_profiles_scan[i] - fixel_mtr_profiles_rescan[i]) / fixel_mtr_profiles_scan[i] * 100

     data_for_boxplot = []
     for sec in range(args.nb_sections):
          # MTR diff values at this section across subjects
          data_for_boxplot.append(mtr_profile_diff[:, sec])
          # Fixel-MTR diff
          data_for_boxplot.append(fixel_mtr_profile_diff[:, sec])

     # Plot profiles
     colors = ['#00A759', '#B45E2F']
     norm = mpl.colors.Normalize(vmin=0.3, vmax=0.7)

     fig = plt.figure(figsize=(8, 7))
     gs = GridSpec(2, 2, width_ratios=[20, 1], height_ratios=[3, 1],
                   wspace=0.15, hspace=0.25)
     
     ax1 = fig.add_subplot(gs[0, 0])

     ax1.fill_between(labels[min_nb_subjects_mask],
                      mtr_profile[min_nb_subjects_mask] - mtr_profile_std[min_nb_subjects_mask],
                      mtr_profile[min_nb_subjects_mask] + mtr_profile_std[min_nb_subjects_mask],
                      color=colors[0], alpha=0.2)
     ax1.fill_between(labels[min_nb_subjects_mask],
                      fixel_mtr_profile[min_nb_subjects_mask] - fixel_mtr_profile_std[min_nb_subjects_mask],
                      fixel_mtr_profile[min_nb_subjects_mask] + fixel_mtr_profile_std[min_nb_subjects_mask],
                      color=colors[1], alpha=0.2)

     ax1.plot(labels[min_nb_subjects_mask], mtr_profile[min_nb_subjects_mask],
              label='MTR', marker='o', color=colors[0])
     ax1.plot(labels[min_nb_subjects_mask],
              fixel_mtr_profile[min_nb_subjects_mask], label='Fixel-wise MTR',
              marker='o', color=colors[1])

     # Add secondary axis for NuFO
     ax2 = ax1.twinx()
     ax2.plot(labels[min_nb_subjects_mask],
              nufo_profile[min_nb_subjects_mask],
              linestyle='--', color="darkgrey", zorder=3)
     ax2.scatter(labels[min_nb_subjects_mask],
                 nufo_profile[min_nb_subjects_mask],
                 c=afd_profile[min_nb_subjects_mask],
                 cmap='Greys', norm=norm, edgecolors='darkgrey', zorder=4)
     complexity_handle = Line2D([0], [0],
                                color='darkgrey', linestyle='--',
                                marker='o', markerfacecolor='white',
                                markeredgecolor='darkgrey',
                                label='Complexity')
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

     # Axis labels and legend
     ax1.set_xlabel('Bundle section')
     ax1.set_ylabel('Mean MTR')
     ax1.set_title('Track-profile of MTR and fixel-wise MTR for the {} bundle'.format(args.in_bundle_name))

     # Combine legends from both axes
     lines_1, labels_1 = ax1.get_legend_handles_labels()
     lines_2, labels_2 = ax2.get_legend_handles_labels()
     ax1.legend(lines_1 + lines_2 + [complexity_handle], labels_1 + labels_2 + ["Complexity"], loc='best')

     ax1.set_ylim(0.33, 0.48)
     ax1.set_xlim(0, args.nb_sections + 1)
     ax1.set_xticks(np.arange(1, args.nb_sections + 1, 1))
     ax2.set_ylim(1, 5)
     ax2.set_yticks(np.arange(1, 6, 1))

     # Absolute difference subplot
     ax3 = fig.add_subplot(gs[1, 0])

     # Create boxplot
     bp = ax3.boxplot(data_for_boxplot, patch_artist=True, showfliers=False,
                      medianprops=dict(color='black'))
     # Color the boxes alternating (MTR = cmap[0], Fixel = cmap[1])
     for i, box in enumerate(bp['boxes']):
          col = colors[0] if i % 2 == 0 else colors[1]
          box.set_facecolor(col)
          # box.set_alpha(0.5)

     positions = np.arange(1.5, 2 * args.nb_sections + 1, 2)
     ax3.set_xticks(positions)
     ax3.set_xticklabels(np.arange(1, args.nb_sections + 1))

     ax3.set_ylabel('|%diff| scan-rescan')
     ax3.set_xlabel('Bundle section')
     ax3.set_ylim(0, 10)
     ax3.set_yticks(np.arange(1, 11, 1))
     ax3.set_xlim(0, 2 * args.nb_sections + 2)
     # ax3.set_xticks(np.arange(1, args.nb_sections + 1, 1))
     # ax3.legend(loc='upper right')
     # ax3.grid(alpha=0.3)

     plt.tight_layout()
     # plt.show()
     plt.savefig(args.out_dir + '/track_profile_{}.png'.format(args.in_bundle_name),
                 dpi=300)


if __name__ == "__main__":
    main()
