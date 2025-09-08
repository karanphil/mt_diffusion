#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import argparse
import logging

from dipy.io.streamline import save_tractogram
import json
import nibabel as nib
import numpy as np

from scilpy.io.streamlines import load_tractogram_with_reference
from scilpy.io.utils import (add_overwrite_arg, add_verbose_arg,
                             assert_inputs_exist, assert_outputs_exist)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('in_tractogram')

    p.add_argument('in_json')

    p.add_argument('out_directory')

    p.add_argument('--in_bundles', nargs='+', required=True,
                   help='List of paths of the bundles (.trk) to analyze.')

    add_verbose_arg(p)
    add_overwrite_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, [args.in_tractogram, args.in_bundles,
                                 args.in_json])
    assert_outputs_exist(parser, args, [args.out_directory])

    logging.info('Loading tractogram...')
    sft = load_tractogram_with_reference(parser, args, args.in_tractogram)
    sft.to_vox()
    sft.to_corner()

    with open(args.in_json, 'r') as f:
        dps_data = json.load(f)

    logging.info('Loading bundles...')
    for b in args.in_bundles:
        bundle = load_tractogram_with_reference(parser, args, b)
        bundle.to_vox()
        bundle.to_corner()
        bundle_name = b.split('/')[-1].split('.trk')[0]
        print(bundle_name)
        indices = np.array(dps_data[bundle_name]['indices'])
        new_bundle = sft[indices]
        save_tractogram(new_bundle, 
                        args.out_directory + '/' + bundle_name + '.trk')


if __name__ == "__main__":
    main()
