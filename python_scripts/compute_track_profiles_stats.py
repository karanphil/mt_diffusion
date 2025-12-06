#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to plot the track-profile of a bundle for both the MTR and the
fixel-wise MTR.
"""

import argparse
import logging
import warnings

import numpy as np
import scipy.stats as stats
from statsmodels.stats import multitest


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

    p.add_argument('--in_mtr_profiles', nargs='+', required=True,
                help='Input bundle-specific MTR track-profile txt files '
                        'for all subjects (scan).')

    p.add_argument('--in_fixel_mtr_profiles', nargs='+', required=True,
                help='Input bundle-specific fixel-wise MTR track-profile '
                        'txt files for all subjects (scan).')

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
                                args.in_fixel_mtr_profiles])
    assert_outputs_exist(parser, args, [args.out_dir])

    # Load profiles
    mtr_profiles= np.array([np.loadtxt(f) for f in args.in_mtr_profiles])
    fixel_mtr_profiles = np.array([np.loadtxt(f) for f in args.in_fixel_mtr_profiles])

    mtr_profiles = np.where(mtr_profiles > 0, mtr_profiles, np.nan)
    fixel_mtr_profiles = np.where(fixel_mtr_profiles > 0,
                                  fixel_mtr_profiles, np.nan)

    # shapiro = stats.shapiro(mtr_profiles, axis=0, nan_policy='omit')
    # for i, pvalue in enumerate(shapiro.pvalue):
    #     if pvalue < 0.05:
    #         print("Section {} of the MTR profiles does not follow a normal distribution. Pvalue = {}".format(i+1, pvalue))

    # shapiro = stats.shapiro(fixel_mtr_profiles, axis=0, nan_policy='omit')
    # for i, pvalue in enumerate(shapiro.pvalue):
    #     if pvalue < 0.05:
    #         print("Section {} of the fixel-MTR profiles does not follow a normal distribution. Pvalue = {}".format(i+1, pvalue))

    # ttest = stats.ttest_ind(mtr_profiles, fixel_mtr_profiles, axis=0, nan_policy='omit', equal_var=False)
    ttest = stats.ttest_rel(mtr_profiles, fixel_mtr_profiles, axis=0, nan_policy='omit')
    # print(ttest.pvalue)
    # _, pvalues_corrected = multitest.fdrcorrection(ttest.pvalue, alpha=0.05, method='indep')
    pvalues_corrected = stats.false_discovery_control(ttest.pvalue, axis=0, method='bh')

    print(pvalues_corrected < 0.05)


if __name__ == "__main__":
    main()
