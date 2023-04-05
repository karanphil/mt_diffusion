import numpy as np
import nibabel as nib
from pathlib import Path
import sys

# À rouler dans le output directory!!!
argv = sys.argv

e1_nifti = Path(argv[1])
wm_mask_nifti = Path(argv[2])

# Load the data
wm_mask_img = nib.load(str(wm_mask_nifti))
e1_img = nib.load(str(e1_nifti))

wm_mask_data = wm_mask_img.get_fdata()
e1_data = e1_img.get_fdata()

sub_ses_dir = Path(e1_nifti.parent.name)
if not sub_ses_dir.is_dir():
    sub_ses_dir.mkdir()

# Define the bins
bin_width = 10
bins = np.arange(0, 90 + bin_width, bin_width)

# Calculate the angle between e1 and B0 field
b0_field = np.array([0, 0, 1])
cos_theta = np.dot(e1_data[..., :3], b0_field)
theta = np.arccos(cos_theta) * 180 / np.pi

wm_mask_bool = (wm_mask_data > 0.9)
for i in range(len(bins) - 1):
    angle_mask_0_90 = (theta >= bins[i]) & (theta < bins[i+1]) 
    angle_mask_90_180 = (180 - theta >= bins[i]) & (180 - theta < bins[i+1])
    angle_mask = angle_mask_0_90 | angle_mask_90_180
    mask = wm_mask_bool & angle_mask
    mask_name = "mask_" + str(bins[i]) + "_to_" + str(bins[i+1]) + "_degrees_angle_bin.nii.gz"
    nib.save(nib.Nifti1Image(mask.astype(np.uint8), e1_img.affine), str(sub_ses_dir / mask_name))