#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to find the labels of crossing bundles based on the number of voxels
share by each section of each bundle.
"""

import argparse
import logging
import json

import nibabel as nib
import numpy as np

from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('--in_labels', nargs='+', required=True,
                   help='Input bundles labels files (nifti). '
                        'One file per bundle.')
    
    p.add_argument('--in_bundles_names', nargs='+',
                   help='Names of the bundles being processed.')
    
    p.add_argument('--min_fraction', type=int, default=0.5,
                   help='Minimum fraction of voxels shared to consider a '
                        'crossing. The number of voxels shared by two '
                        'section of two bundles is compared to the total '
                        'number of voxels in each section. Default is 0.5.')

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

    assert_inputs_exist(parser, [args.in_labels])

    if args.in_bundles_names:
        if len(args.in_bundles_names) != len(args.in_labels):
            parser.error('Number of bundle names should match number of '
                         'labels files.')
        bundle_names = args.in_bundles_names
    else:
        # extract names from filenames (labels_bundlename.nii.gz)
        bundle_names = []
        for label_file in args.in_labels:
            base_name = label_file.split('/')[-1]
            bundle_names.append(base_name.replace('labels_', '').replace('.nii.gz', ''))

    # Load labels
    labels_data = []
    for label_file in args.in_labels:
        img = nib.load(label_file)
        labels_data.append(img.get_fdata().astype(np.int16))
    
    nb_bundles = len(labels_data)
    nb_sections = args.nb_sections
    
    # Compute crossing bundles labels
    crossing_info = {}
    for i in range(nb_bundles):
        bundle_i_name = bundle_names[i]
        crossing_info[bundle_i_name] = {}
        print(f'Processing bundle {bundle_i_name}')
        for section_i in range(1, nb_sections + 1):
            section_i_mask = (labels_data[i] == section_i)
            crossing_info[bundle_i_name][section_i] = [("Nb_voxels", np.sum(section_i_mask))]
            for j in range(nb_bundles):
                if i == j:
                    continue
                bundle_j_name = bundle_names[j]
                # print(f'  Comparing with bundle {bundle_j_name}')
                for section_j in range(1, nb_sections + 1):
                    section_j_mask = (labels_data[j] == section_j)
                    shared_voxels = np.sum(section_i_mask & section_j_mask)
                    # if (shared_voxels >= args.min_fraction * np.sum(section_i_mask)) \
                        # or (shared_voxels >= args.min_fraction * np.sum(section_j_mask)):
                    crossing_info[bundle_i_name][section_i].append(
                        (bundle_j_name, section_j, shared_voxels))
                        
    # # Print crossing bundles labels
    # for bundle_name, sections in crossing_info.items():
    #     print(f'Crossing bundles for {bundle_name}:')
    #     for section, crossings in sections.items():
    #         if crossings:
    #             crossing_str = ', '.join(
    #                 [f'{bname} (section {s}, voxels {nvox})' 
    #                  for bname, s, nvox in crossings])
    #             print(f'  Section {section}: {crossing_str}')
    #     print('')

    # Saving results in json file
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super(NpEncoder, self).default(obj)
    out_json = 'crossing_bundles_labels.json'
    with open(out_json, 'w') as f:
        json.dump(crossing_info, f, indent=4, cls=NpEncoder)


if __name__ == "__main__":
    main()
