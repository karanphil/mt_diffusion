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

# bins_width = [1, 3]
bins_width = [10]
norm_thr = 0.7

for w in bins_width: # width of the angle bins

    # Define the bins
    # bins = np.arange(0, 180, 5)
    bins = np.arange(0, 90 + w, w)
    mid_bins = (bins[:-1] + bins[1:]) / 2.

    # Calculate the angle between e1 and B0 field
    b0_field = np.array([0, 0, 1])
    cos_theta_f1 = np.dot(peaks_data[..., 0:3], b0_field)
    theta_f1 = np.arccos(cos_theta_f1) * 180 / np.pi
    cos_theta_f2 = np.dot(peaks_data[..., 3:6], b0_field)
    theta_f2 = np.arccos(cos_theta_f2) * 180 / np.pi

    mtr_means = np.zeros((len(bins) - 1, len(bins) - 1))
    ihmtr_means = np.zeros((len(bins) - 1, len(bins) - 1))
    nb_voxels = np.zeros((len(bins) - 1, len(bins) - 1))

    # Apply the WM mask
    wm_mask_bool = (wm_mask_data > 0.9) & (nufo_data == 2)
    norm_mask_bool = (peaks_norm[..., 0] > norm_thr) & (peaks_norm[..., 1] > norm_thr)

    for i in range(len(bins) - 1):
        angle_mask_0_90 = (theta_f1 >= bins[i]) & (theta_f1 < bins[i+1])
        angle_mask_90_180 = (180 - theta_f1 >= bins[i]) & (180 - theta_f1 < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask_f1 = angle_mask
        for j in range(len(bins) - 1):
            angle_mask_0_90 = (theta_f2 >= bins[j]) & (theta_f2 < bins[j+1]) 
            angle_mask_90_180 = (180 - theta_f2 >= bins[j]) & (180 - theta_f2 < bins[j+1])
            angle_mask = angle_mask_0_90 | angle_mask_90_180
            mask_f2 = angle_mask
            mask = mask_f1 & mask_f2 & wm_mask_bool & norm_mask_bool
            nb_voxels[i, j] = np.sum(mask)
            if np.sum(mask) < 5:
                mtr_means[i, j] = None
                ihmtr_means[i, j] = None
            else:
                mtr_means[i, j] = np.mean(mtr_data[mask])
                ihmtr_means[i, j] = np.mean(ihmtr_data[mask])
    
    for i in range(len(bins) - 1):
        for j in range(i):
            mtr_means[i, j] = (mtr_means[i, j] + mtr_means[j, i]) / 2
            mtr_means[j, i] = mtr_means[i, j]
            ihmtr_means[i, j] = (ihmtr_means[i, j] + ihmtr_means[j, i]) / 2
            ihmtr_means[j, i] = ihmtr_means[i, j]
            nb_voxels[i, j] = nb_voxels[i, j] + nb_voxels[j, i]
            nb_voxels[j, i] = nb_voxels[i, j]

    mtr_means_diag = np.diagonal(mtr_means)
    ihmtr_means_diag = np.diagonal(ihmtr_means)
    nb_voxels_diag = np.diagonal(nb_voxels)

    # Save the results to a text file
    results = np.column_stack((bins[:-1], bins[1:], mtr_means, ihmtr_means, nb_voxels))
    txt_name = "results_" + str(w) + "_degrees_bins.txt"
    txt_path = sub_ses_dir / txt_name
    # !!!!!!!! À sauver en disant 3D ou 2 fibers ou ailleurs!!!!!!!
    # if ratios:
    #     np.savetxt(str(txt_path), results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTR_mean\tihMTR_mean\tNb_voxels')
    # elif sats:
    #     np.savetxt(str(txt_path), results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTsat_mean\tihMTsat_mean\tNb_voxels')

    # Plot diagonal
    max_count = np.max(nb_voxels_diag)
    norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
    plot_init()
    fig, (ax1, ax2, cax) = plt.subplots(1, 3, gridspec_kw={"width_ratios":[1, 1, 0.05]})
    ax1.scatter(mid_bins, mtr_means_diag, c=nb_voxels_diag, cmap='Greys', norm=norm, edgecolors='black', linewidths=1)
    colorbar = ax1.set_xlabel('Angle between e1 and B0 field (degrees)')
    if ratios:
        ax1.set_ylabel('MTR mean')
        ax1.set_title('MTR vs Angle')
    elif sats:
        ax1.set_ylabel('MTsat mean')
        ax1.set_title('MTsat vs Angle')
    colorbar = ax2.scatter(mid_bins, ihmtr_means_diag, c=nb_voxels_diag, cmap='Greys', norm=norm, edgecolors='black', linewidths=1)
    ax2.set_xlabel('Angle between e1 and B0 field (degrees)')
    if ratios:
        ax2.set_ylabel('ihMTR mean')
        ax2.set_title('ihMTR vs Angle')
    elif sats:
        ax2.set_ylabel('ihMTsat mean')
        ax2.set_title('ihMTsat vs Angle')
    fig.colorbar(colorbar, cax=cax, label="Voxel count")
    fig.tight_layout()
    plt.show()
    fig_name = "plot_" + str(j) + "_degrees_bins.png"
    fig_path = sub_ses_dir / fig_name
    # plt.savefig(fig_path, dpi=300)
    plt.close()

    # Plot fiber 2 for fiber 1 fixed.
    plot_init()
    fig, (ax1, ax2) = plt.subplots(1, 2)
    for idx in range(mid_bins.shape[0]):
        ax1.plot(mid_bins, mtr_means[idx, :], "o-", color="C" + str(idx), label=str(mid_bins[idx]) + " degrees")
        ax2.plot(mid_bins, ihmtr_means[idx, :], "o-", color="C" + str(idx), label=str(mid_bins[idx]) + " degrees")

    ax1.set_xlabel('Angle of fiber 2')
    if ratios:
        ax1.set_ylabel('MTR mean')
        ax1.set_title('MTR vs Angle')
    elif sats:
        ax1.set_ylabel('MTsat mean')
        ax1.set_title('MTsat vs Angle')
    ax2.set_xlabel('Angle of fiber 2')
    if ratios:
        ax2.set_ylabel('ihMTR mean')
        ax2.set_title('ihMTR vs Angle')
    elif sats:
        ax2.set_ylabel('ihMTsat mean')
        ax2.set_title('ihMTsat vs Angle')
    ax2.legend()
    fig.tight_layout()
    plt.show()
    fig_name = "plot_" + str(j) + "_degrees_bins.png"
    fig_path = sub_ses_dir / fig_name
    # plt.savefig(fig_path, dpi=300)
    plt.close()

    # Plot the results
    print("Step 3:")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(mid_bins, mid_bins)
    ax.plot_surface(X, Y, mtr_means)
    ax.set_xlabel('Angle of fiber 1')
    ax.set_ylabel('Angle of fiber 2')
    ax.set_zlabel('MTR mean')
    ax.set_title('MTR vs Angle')
    fig.tight_layout()
    plt.show()
    plt.close()