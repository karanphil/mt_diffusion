import numpy as np
import nibabel as nib
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# À rouler dans le output directory!!!
argv = sys.argv

# Select MTR and ihMTR or MTsat and ihMTsat
# ratios = False
# sats = True
ratios = True
sats = False

def plot_init():
    # plt.rcParams["font.family"] = "serif"
    # plt.rcParams['font.serif'] = 'Helvetica'
    # plt.style.use('seaborn-notebook')
    plt.rcParams['axes.grid'] = False
    plt.rcParams['grid.color'] = "darkgrey"
    plt.rcParams['grid.linewidth'] = 1
    plt.rcParams['grid.linestyle'] = "-"
    plt.rcParams['grid.alpha'] = "0.5"
    plt.rcParams['figure.figsize'] = (10.0, 5.0)
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.2*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.linewidth'] =1
    plt.rcParams['lines.linewidth']=1
    plt.rcParams['lines.markersize']=4
    # plt.rcParams['text.latex.unicode']=True
    # plt.rcParams['text.latex.preamble'] = [r'\usepackage{amssymb}', r"\usepackage{amstext}"]
    # plt.rcParams['mathtext.default']='regular'

peaks_nifti = Path(argv[1])
mtr_nifti = Path(argv[2])
ihmtr_nifti = Path(argv[3])
wm_mask_nifti = Path(argv[4])
nufo_nifti = Path(argv[5])

# Load the data
mtr_img = nib.load(str(mtr_nifti))
ihmtr_img = nib.load(str(ihmtr_nifti))
wm_mask_img = nib.load(str(wm_mask_nifti))
peaks_img = nib.load(str(peaks_nifti))
nufo_img = nib.load(str(nufo_nifti))

mtr_data = mtr_img.get_fdata()
ihmtr_data = ihmtr_img.get_fdata()
wm_mask_data = wm_mask_img.get_fdata()
peaks_data = peaks_img.get_fdata()
nufo_data = nufo_img.get_fdata()

sub_ses_dir = Path(peaks_nifti.parent.name)
if not sub_ses_dir.is_dir():
    sub_ses_dir.mkdir()

nb_peaks = int(peaks_data.shape[-1] / 3)
peaks_norm = np.zeros((peaks_data.shape[:-1] + (nb_peaks,)))
for i in range(nb_peaks):
    peaks_norm[..., i] = np.linalg.norm(peaks_data[..., i * 3 : (i + 1) * 3], axis=-1)
    for j in range(3):
        peaks_data[..., i * 3 + j] /= peaks_norm[..., i]

bins_width = [1, 3]
# bins_width = [1]

for j in bins_width: # width of the angle bins

    # Define the bins
    # bins = np.arange(0, 180, 5)
    bins = np.arange(0, 90 + j, j)

    # Calculate the angle between e1 and B0 field
    b0_field = np.array([0, 0, 1])
    cos_theta_f1 = np.dot(peaks_data[..., 0:3], b0_field)
    theta_f1 = np.arccos(cos_theta_f1) * 180 / np.pi
    cos_theta_f2 = np.dot(peaks_data[..., 3:6], b0_field)
    theta_f2 = np.arccos(cos_theta_f2) * 180 / np.pi

    mtr_means = np.zeros((2, len(bins) - 1))
    ihmtr_means = np.zeros((2, len(bins) - 1))
    nb_voxels = np.zeros((2, len(bins) - 1))

    # Apply the WM mask
    wm_mask_bool = (wm_mask_data > 0.9) &  (nufo_data == 2)
    norm_mask_bool = (peaks_norm[..., 0] > 0.9) & (peaks_norm[..., 1] > 0.9)

    for f in range(2):
        if f == 0:
            theta = theta_f1
        elif f == 1:
            theta = theta_f2
        for i in range(len(bins) - 1):
            angle_mask_0_90 = (theta >= bins[i]) & (theta < bins[i+1]) 
            angle_mask_90_180 = (180 - theta >= bins[i]) & (180 - theta < bins[i+1])
            angle_mask = angle_mask_0_90 | angle_mask_90_180
            mask = wm_mask_bool & angle_mask
            nb_voxels[f, i] = np.sum(mask)
            if np.sum(mask) < 5:
                mtr_means[f, i] = None
                ihmtr_means[f, i] = None
            else:
                mtr_means[f, i] = np.mean(mtr_data[mask])
                ihmtr_means[f, i] = np.mean(ihmtr_data[mask])

        # Save the results to a text file
        results = np.column_stack((bins[:-1], bins[1:], mtr_means[k], ihmtr_means[k], nb_voxels[k]))
        txt_name = "results_" + str(j) + "_degrees_bins.txt"
        txt_path = sub_ses_dir / txt_name
        if ratios:
            np.savetxt(str(txt_path), results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTR_mean\tihMTR_mean\tNb_voxels')
        elif sats:
            np.savetxt(str(txt_path), results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTsat_mean\tihMTsat_mean\tNb_voxels')

    # Plot the results
    max_count = np.max(nb_voxels[:])
    norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
    plot_init()
    fig, (ax1, ax2, cax) = plt.subplots(1, 3, gridspec_kw={"width_ratios":[1,1, 0.05]})
    ax1.scatter(bins[:-1], mtr_means[:], c=nb_voxels[:], cmap='Greys', norm=norm, edgecolors='black', linewidths=1)
    ax1.set_xlabel('Angle between e1 and B0 field (degrees)')
    if ratios:
        ax1.set_ylabel('MTR mean')
        ax1.set_title('MTR vs Angle')
    elif sats:
        ax1.set_ylabel('MTsat mean')
        ax1.set_title('MTsat vs Angle')
    colorbar = ax2.scatter(bins[:-1], ihmtr_means[:], c=nb_voxels[:], cmap='Greys', norm=norm, edgecolors='black', linewidths=1)
    ax2.set_xlabel('Angle between e1 and B0 field (degrees)')
    if ratios:
        ax2.set_ylabel('ihMTR mean')
        ax2.set_title('ihMTR vs Angle')
    elif sats:
        ax2.set_ylabel('ihMTsat mean')
        ax2.set_title('ihMTsat vs Angle')
    fig.colorbar(colorbar, cax=cax, label="Voxel count")
    fig.tight_layout()
    # plt.show()
    fig_name = "plot_" + str(j) + "_degrees_bins.png"
    fig_path = sub_ses_dir / fig_name
    plt.savefig(fig_path, dpi=300)
    plt.close()