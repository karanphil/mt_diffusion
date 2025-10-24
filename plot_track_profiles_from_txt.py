#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to plot the track-profile of a bundle for both the MTR and the
fixel-wise MTR.
"""

import argparse
import logging
import warnings

from cmcrameri import cm
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

    p.add_argument('out_dir')

    p.add_argument('--in_mtr_profiles', nargs='+', required=True,
                   help='Input MTR profile text files.')
    
    p.add_argument('--in_fixel_mtr_profiles', nargs='+', required=True,
                   help='Input fixel-wise MTR profile text files.')
    
    p.add_argument('--in_nufo_profiles', nargs='+', required=True,
                   help='Input NuFO profile text files.')
    
    p.add_argument('--in_afd_profiles', nargs='+', required=True,
                   help='Input AFD profile text files.')

    p.add_argument('--bundle_name')

    p.add_argument('--variance', action='store_true',
                   help='Compute and display the variance between profiles of '
                        'the MTR and fixel-wise MTR for each bundle section. '
                        'Default is False.')
    
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

    assert_inputs_exist(parser, [args.in_mtr_profiles,
                                 args.in_fixel_mtr_profiles,
                                 args.in_nufo_profiles, args.in_afd_profiles])
    assert_outputs_exist(parser, args, [args.out_dir])

    # Load profiles
    mtr_profiles = np.array([np.loadtxt(f) for f in args.in_mtr_profiles])
    fixel_mtr_profiles = np.array([np.loadtxt(f) for f in args.in_fixel_mtr_profiles])
    nufo_profiles = np.array([np.loadtxt(f) for f in args.in_nufo_profiles])
    afd_profiles = np.array([np.loadtxt(f) for f in args.in_afd_profiles])
    labels = np.arange(1, len(mtr_profiles[0]) + 1, 1)

    mtr_profiles = np.where(mtr_profiles > 0, mtr_profiles, np.nan)
    fixel_mtr_profiles = np.where(fixel_mtr_profiles > 0, fixel_mtr_profiles,
                                  np.nan)
    nufo_profiles = np.where(nufo_profiles > 0, nufo_profiles, np.nan)
    afd_profiles = np.where(afd_profiles > 0, afd_profiles, np.nan)

    nb_subjects = np.sum(~np.isnan(fixel_mtr_profiles), axis=0)
    min_nb_subjects_mask = nb_subjects >= args.min_nb_subjects

    mtr_profile = np.nanmean(mtr_profiles, axis=0)
    fixel_mtr_profile = np.nanmean(fixel_mtr_profiles, axis=0)
    nufo_profile = np.nanmean(nufo_profiles, axis=0)
    afd_profile = np.nanmean(afd_profiles, axis=0)
    mtr_profile_std = np.nanstd(mtr_profiles, axis=0)
    fixel_mtr_profile_std = np.nanstd(fixel_mtr_profiles, axis=0)

    mtr_profile = np.where(np.isnan(mtr_profile), 0, mtr_profile)
    fixel_mtr_profile = np.where(np.isnan(fixel_mtr_profile), 0, fixel_mtr_profile)
    nufo_profile = np.where(np.isnan(nufo_profile), 0, nufo_profile)
    afd_profile = np.where(np.isnan(afd_profile), 0, afd_profile)

    # Plot profiles
    cmap = cm.naviaS
    cmap_idx = np.arange(3, 1000, 1)
    norm = mpl.colors.Normalize(vmin=0.3, vmax=0.7)

    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(labels[min_nb_subjects_mask], mtr_profile[min_nb_subjects_mask],
             label='MTR', marker='o', color=cmap(cmap_idx[0]))
    ax1.plot(labels[min_nb_subjects_mask],
             fixel_mtr_profile[min_nb_subjects_mask], label='Fixel-wise MTR',
             marker='o', color=cmap(cmap_idx[1]))

    # import matplotlib.colors as mcolors
    # colors = mcolors.TABLEAU_COLORS
    # for i in range(fixel_mtr_profiles.shape[0]):
    #     ax1.scatter(labels[fixel_mtr_profiles[i] != 0],
    #                 fixel_mtr_profiles[i][fixel_mtr_profiles[i] != 0],
    #                 zorder=1)

    if args.variance:
        ax1.fill_between(labels[min_nb_subjects_mask],
                            mtr_profile[min_nb_subjects_mask] - mtr_profile_std[min_nb_subjects_mask],
                            mtr_profile[min_nb_subjects_mask] + mtr_profile_std[min_nb_subjects_mask],
                            color=cmap(cmap_idx[0]), alpha=0.2)
        ax1.fill_between(labels[min_nb_subjects_mask],
                            fixel_mtr_profile[min_nb_subjects_mask] - fixel_mtr_profile_std[min_nb_subjects_mask],
                            fixel_mtr_profile[min_nb_subjects_mask] + fixel_mtr_profile_std[min_nb_subjects_mask],
                            color=cmap(cmap_idx[1]), alpha=0.2)

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
    sm = plt.cm.ScalarMappable(cmap='Greys',
                               norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax2, pad=0.08)
    cbar.set_label('Mean AFD', color="darkgrey")
    cbar.ax.tick_params(color='darkgrey', labelcolor='darkgrey')
    cbar.outline.set_edgecolor('darkgrey')
    # cbar.tick_params(axis='y', labelcolor="darkgrey")

    # Axis labels and legend
    ax1.set_xlabel('Bundle section')
    ax1.set_ylabel('Mean MTR')
    bundle_name = f' for the {args.bundle_name} bundle' if args.bundle_name else ''
    ax1.set_title(f'Track-profile of MTR and fixel-wise MTR {bundle_name}')

    # Combine legends from both axes
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2 + [complexity_handle], labels_1 + labels_2 + ["Complexity"], loc='best')

    ax1.set_ylim(0.33, 0.45)
    ax1.set_xlim(0, 21)
    ax1.set_xticks(np.arange(1, 21, 1))
    ax2.set_ylim(1, 5)
    ax2.set_yticks(np.arange(1, 6, 1))

    plt.tight_layout()
    # plt.show()
    plt.savefig(args.out_dir + '/track_profile_{}.png'.format(args.bundle_name),
                dpi=300)


if __name__ == "__main__":
    main()
