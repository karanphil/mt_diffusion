import argparse

import nibabel as nib
import numpy as np


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('in_voxel_density_masks')

    p.add_argument('in_LUT')

    p.add_argument('out_dir')

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    voxel_density_img = nib.load(args.in_voxel_density_masks)
    bundles_mask = voxel_density_img.get_fdata().astype(bool)

    lut = np.loadtxt(args.in_LUT, dtype=str)

    for i in range(bundles_mask.shape[-1]):

        bundle_name = lut[0][i]
        bundle_mask = bundles_mask[..., i].astype(np.uint8)

        nib.save(nib.Nifti1Image(bundle_mask, voxel_density_img.affine), args.out_dir + f'/voxel_density_mask_voxel-norm_{bundle_name}.nii.gz')


if __name__ == "__main__":
    main()
