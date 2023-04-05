import numpy as np
import nibabel as nib
from pathlib import Path
import sys

# À rouler dans le output directory!!!
argv = sys.argv

wm_mask_nifti = Path(argv[5])
e1_nifti = Path(argv[1])

# Load the data
wm_mask_img = nib.load(str(wm_mask_nifti))
e1_img = nib.load(str(e1_nifti))

wm_mask_data = wm_mask_img.get_fdata()
e1_data = e1_img.get_fdata()

sub_ses_dir = Path(e1_nifti.parent.name)
if not sub_ses_dir.is_dir():
    sub_ses_dir.mkdir()

# Define the bins
bin_width = 20
bins = np.arange(0, 180 + bin_width, bin_width)

# Calculate the angle between e1 and B0 field
b0_field = np.array([0, 0, 1])
cos_theta = np.dot(e1_data[..., :3], b0_field)
theta = np.arccos(cos_theta) * 180 / np.pi

wm_mask_bool = (wm_mask_data > 0)
for i in range(len(bins) - 1):
    angle_mask = (theta >= bins[i]) & (theta < bins[i+1])
    mask = wm_mask_bool & angle_mask