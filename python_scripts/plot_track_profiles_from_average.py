#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to plot the track-profile of a bundle for both the MTR and the
fixel-wise MTR. Made for a single profile at a time.
"""

import argparse
import logging

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

    p.add_argument('in_mtr_profile',
                   help='Input bundle-specific MTR track-profile txt file.')

    p.add_argument('in_fixel_mtr_profile',
                   help='Input bundle-specific fixel-wise MTR track-profile '
                        'txt file.')

    p.add_argument('in_bundle_name', type=str,
                   help='Name of the bundle being processed.')

    p.add_argument('in_afd_fixel_profile',
                   help='Input bundle-specific AFD fixel track-profile txt '
                        'file.')

    p.add_argument('in_nufo_profile',
                   help='Input NuFO track-profile txt file.')

    p.add_argument('out_dir',
                   help='Output directory to save the track-profiles plot.')
    
    p.add_argument('--nb_sections', type=int, default=20,
                   help='Number of sections in the bundle labels. '
                        'Default is 20.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_mtr_profile,
                                 args.in_fixel_mtr_profile,
                                 args.in_afd_fixel_profile,
                                 args.in_nufo_profile])
    assert_outputs_exist(parser, args, [args.out_dir])

    # Load profiles
    mtr_profile = np.array(np.loadtxt(args.in_mtr_profile))
    fixel_mtr_profile = np.array(np.loadtxt(args.in_fixel_mtr_profile))
    afd_profile = np.array(np.loadtxt(args.in_afd_fixel_profile))
    nufo_profile = np.array(np.loadtxt(args.in_nufo_profile))
    labels = np.arange(1, args.nb_sections + 1, 1)

    print(mtr_profile.shape)

    colors = ['#00A759', '#B45E2F']
    norm = mpl.colors.Normalize(vmin=0.3, vmax=0.7)

    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(labels[mtr_profile != 0], mtr_profile[mtr_profile != 0],
             label='MTR', marker='o', color=colors[0])
    ax1.plot(labels[fixel_mtr_profile != 0],
             fixel_mtr_profile[fixel_mtr_profile != 0], label='Fixel-wise MTR',
             marker='o', color=colors[1])

    # Add secondary axis for NuFO
    ax2 = ax1.twinx()
    ax2.plot(labels[nufo_profile != 0],
             nufo_profile[nufo_profile != 0],
             linestyle='--', color="darkgrey", zorder=3)
    ax2.scatter(labels[nufo_profile != 0],
                nufo_profile[nufo_profile != 0],
                c=afd_profile[nufo_profile != 0],
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
    ax1.set_title('Track-profile of MTR and fixel-wise MTR for the {} bundle'.format(args.in_bundle_name))

    # Combine legends from both axes
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2 + [complexity_handle], labels_1 + labels_2 + ["Complexity"], loc='best')

    ax1.set_ylim(0.33, 0.45)
    ax1.set_xlim(0, args.nb_sections + 1)
    ax1.set_xticks(np.arange(1, args.nb_sections + 1, 1))
    ax2.set_ylim(1, 5)
    ax2.set_yticks(np.arange(1, 6, 1))

    plt.tight_layout()
    # plt.show()
    plt.savefig(args.out_dir + '/track_profile_{}.png'.format(args.in_bundle_name),
                dpi=300)


if __name__ == "__main__":
    main()
