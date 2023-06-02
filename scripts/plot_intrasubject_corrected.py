import argparse
import numpy as np
import nibabel as nib
from pathlib import Path

from modules.io import plot_multiple_means

from scilpy.io.utils import (add_overwrite_arg)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('out_name',
                   help='Path of the output name.')
    
    p.add_argument('--in_results', nargs='+',
                   help='Path to the results txt.')
    p.add_argument('--in_cr_results', nargs='+',
                   help='Path to the corrected results txt.')
    p.add_argument('--input_type', default="ratios",
                   help='Type of input measures (ratios or sats).')
    add_overwrite_arg(p)
    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    in_results = args.in_results[0:5]

    dims = np.loadtxt(in_results[0], skiprows=1).shape
    mt_means = np.zeros((len(in_results), dims[0]))
    ihmt_means = np.zeros((len(in_results), dims[0]))
    nb_voxels = np.zeros((len(in_results), dims[0]))
    labels = np.zeros((len(in_results)), dtype=object)

    for i, result in enumerate(in_results):
        results = np.loadtxt(result, skiprows=1)
        mt_means[i] = results[:, 2]
        ihmt_means[i] = results[:, 3]
        nb_voxels[i] = results[:, 4]

        labels[i] = "Session " + str(i+1)

    bin_min = results[:, 0]
    bin_max = results[:, 1]
    bins = np.zeros((bin_min.shape[0] + 1))
    bins[:-1] = bin_min
    bins[-1] = bin_max[-1]

    if args.in_cr_results:
        in_cr_results = args.in_cr_results[0:5]
        mt_cr_means = np.zeros((len(in_cr_results), dims[0]))
        ihmt_cr_means = np.zeros((len(in_cr_results), dims[0]))

        for i, result in enumerate(in_cr_results):
            results = np.loadtxt(result, skiprows=1)
            mt_cr_means[i] = results[:, 2]
            ihmt_cr_means[i] = results[:, 3]
    else:
        mt_cr_means = None
        ihmt_cr_means = None

    plot_multiple_means(bins, mt_means, ihmt_means, nb_voxels, args.out_name,
                        labels, mt_cr_means=mt_cr_means,
                        ihmt_cr_means=ihmt_cr_means, input_dtype=args.input_type)


if __name__ == "__main__":
    main()
