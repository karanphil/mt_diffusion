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
    
    p.add_argument('--in_subjects', nargs='+',
                   help='List of all subjects.')
    p.add_argument('--in_name',
                   help='Name of the results txt.')
    p.add_argument('--in_name_cr',
                   help='Name of the corrected results txt.')
    p.add_argument('--input_type', default="ratios",
                   help='Type of input measures (ratios or sats).')
    add_overwrite_arg(p)
    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    subjects = args.in_subjects

    dummy = subjects[0] + "_ses-1/"
    dims = np.loadtxt(dummy + args.in_name, skiprows=1).shape
    results = np.zeros(((len(subjects),) + dims))
    labels = np.zeros((len(subjects)), dtype=object)

    for i, subject in enumerate(subjects):
        print(subject)
        sessions = list(Path('.').glob(subject + "*"))
        for session in sessions[:5]:
            print(session)
            results[i] += np.loadtxt(str(session) + "/" + args.in_name, skiprows=1)
        results[i] /= 5
        labels[i] = str(i+1)

    bins = np.zeros((results[:, :, 0].shape[1] + 1))
    bins[:results[:, :, 0].shape[1]] = results[0, :, 0]
    bins[-1] = results[0, -1, 1]

    mt_means = results[:, :, 2]
    ihmt_means = results[:, :, 3]
    nb_voxels = results[:, :, 4]

    if args.in_name_cr:
        results_cr = np.zeros(((len(subjects),) + dims))
        print(subject)
        sessions = list(Path('.').glob(subject + "*"))
        for session in sessions[:5]:
            print(session)
            results_cr[i] += np.loadtxt(str(session) + "/" + args.in_name_cr, skiprows=1)
        results_cr[i] /= 5
        mt_cr_means = results_cr[:, :, 2]
        ihmt_cr_means = results_cr[:, :, 3]
    else:
        mt_cr_means = None
        ihmt_cr_means = None

    plot_multiple_means(bins, mt_means, ihmt_means, nb_voxels, args.out_name,
                        mt_cr_means=mt_cr_means,
                        ihmt_cr_means=ihmt_cr_means, input_dtype=args.input_type)


if __name__ == "__main__":
    main()
