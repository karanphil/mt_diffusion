import numpy as np
import nibabel as nib
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# À rouler dans le output directory!!!
argv = sys.argv

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
    plt.rcParams['lines.markersize']=3
    # plt.rcParams['text.latex.unicode']=True
    # plt.rcParams['text.latex.preamble'] = [r'\usepackage{amssymb}', r"\usepackage{amstext}"]
    # plt.rcParams['mathtext.default']='regular'

e1_nifti = Path(argv[1])
fa_nifti = Path(argv[2])
mtr_nifti = Path(argv[3])
ihmtr_nifti = Path(argv[4])
wm_mask_nifti = Path(argv[5])
nufo_nifti = Path(argv[6])

# Load the data
mtr_img = nib.load(str(mtr_nifti))
ihmtr_img = nib.load(str(ihmtr_nifti))
fa_img = nib.load(str(fa_nifti))
wm_mask_img = nib.load(str(wm_mask_nifti))
e1_img = nib.load(str(e1_nifti))
nufo_img = nib.load(str(nufo_nifti))

mtr_data = mtr_img.get_fdata()
ihmtr_data = ihmtr_img.get_fdata()
fa_data = fa_img.get_fdata()
wm_mask_data = wm_mask_img.get_fdata()
e1_data = e1_img.get_fdata()
nufo_data = nufo_img.get_fdata()

sub_ses_dir = Path(fa_nifti.parent.name)
if not sub_ses_dir.is_dir():
    sub_ses_dir.mkdir()

# bins_width = [1, 3]
# fa_thrs = [0.5, 0.6, 0.7]
bins_width = [1]
fa_thrs = [0.5]

for j in bins_width: # width of the angle bins

    # Define the bins
    # bins = np.arange(0, 180, 5)
    bins = np.arange(0, 180 + j, j)

    # Calculate the angle between e1 and B0 field
    b0_field = np.array([0, 0, 1])
    cos_theta = np.dot(e1_data[..., :3], b0_field)
    theta = np.arccos(cos_theta) * 180 / np.pi

    for l in [True, False]: # use NuFo or not
        mtr_means = np.zeros((len(fa_thrs), len(bins) - 1))
        ihmtr_means = np.zeros((len(fa_thrs), len(bins) - 1))
        nb_voxels = np.zeros((len(fa_thrs), len(bins) - 1))

        for k, fa_thr in enumerate(fa_thrs): # FA threshold

            # Apply the WM mask and FA threshold
            if l:
                wm_mask_bool = (wm_mask_data > 0) & (fa_data > fa_thr) & (nufo_data == 1)
            else:
                wm_mask_bool = (wm_mask_data > 0) & (fa_data > fa_thr)

            for i in range(len(bins) - 1):
                angle_mask = (theta >= bins[i]) & (theta < bins[i+1])
                mask = wm_mask_bool & angle_mask
                nb_voxels[k, i] = np.sum(mask)
                if np.sum(mask) < 5:
                    mtr_means[k, i] = None
                    ihmtr_means[k, i] = None
                else:
                    mtr_means[k, i] = np.mean(mtr_data[mask])
                    ihmtr_means[k, i] = np.mean(ihmtr_data[mask])

            # Save the results to a text file
            results = np.column_stack((bins[:-1], bins[1:], mtr_means[k], ihmtr_means[k], nb_voxels[k]))
            txt_name = "results_" + str(j) + "_degrees_bins_" + str(fa_thr) + "_FA_thr_NuFo_" + str(l) + ".txt"
            txt_path = sub_ses_dir / txt_name
            # np.savetxt(str(txt_path), results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTR_mean\tihMTR_mean\tNb_voxels')

            for i in [90, 180]: # range of the angle bins to visualize
                # Plot the results
                if i == 90:
                    angles = int((len(bins) - 1) / 2) + 1
                    means = int((len(bins) - 1) / 2) + 1
                elif i == 180:
                    angles = -1
                    means = mtr_means.shape[-1]
                plot_init()
                fig, (ax1, ax2, cax) = plt.subplots(1, 3, gridspec_kw={"width_ratios":[1,1, 0.05]})
                fig.subplots_adjust(wspace=0.3)
                ax1.scatter(bins[:angles], mtr_means[k, :means], c=nb_voxels[k, :means], cmap='Greens', edgecolors='black', linewidths=0.5)
                ax1.set_xlabel('Angle between e1 and B0 field (degrees)')
                ax1.set_ylabel('MTR mean')
                ax1.set_title('MTR vs Angle')
                # ax3 = ax1.twinx()
                # ax3.scatter(bins[:angles], nb_voxels[k, :means])
                # ax3.set_ylabel('Voxel count')
                # plt.subplot(1, 2, 2)
                colorbar = ax2.scatter(bins[:angles], ihmtr_means[k, :means], c=nb_voxels[k, :means], cmap='Greens', edgecolors='black', linewidths=0.5)
                ax2.set_xlabel('Angle between e1 and B0 field (degrees)')
                ax2.set_ylabel('ihMTR mean')
                ax2.set_title('ihMTR vs Angle')
                # ax4 = ax2.twinx()
                # ax4.scatter(bins[:angles], nb_voxels[k, :means])
                # ax4.set_ylabel('Voxel count')
                fig.colorbar(colorbar, cax=cax, label="Voxel count")
                fig.tight_layout()
                plt.show()
                fig_name = "plot_" + str(j) + "_degrees_bins_" + str(fa_thr) + "_FA_thr_NuFo_" + str(l) + "_" + str(i) + "_degrees_range" + ".png"
                fig_path = sub_ses_dir / fig_name
                #plt.savefig(fig_path, dpi=300)
                plt.close()
                # plt.figure(figsize=(10, 5))
                # plt.subplot(1, 2, 1)
                # plt.scatter(bins[:angles], mtr_means[k, :means], 2)
                # plt.xlabel('Angle between e1 and B0 field (degrees)')
                # plt.ylabel('MTR mean')
                # plt.title('MTR vs Angle')
                # plt.subplot(1, 2, 2)
                # plt.scatter(bins[:angles], ihmtr_means[k, :means], 2)
                # plt.xlabel('Angle between e1 and B0 field (degrees)')
                # plt.ylabel('ihMTR mean')
                # plt.title('ihMTR vs Angle')
                # plt.tight_layout()
                # # plt.show()
                # fig_name = "plot_" + str(j) + "_degrees_bins_" + str(fa_thr) + "_FA_thr_NuFo_" + str(l) + "_" + str(i) + "_degrees_range" + ".png"
                # fig_path = sub_ses_dir / fig_name
                # plt.savefig(fig_path)
                # plt.close()

        # for i in [90, 180]: # range of the angle bins to visualize
        #     # Plot the results
        #     if i == 90:
        #         angles = int((len(bins) - 1) / 2) + 1
        #         means = int((len(bins) - 1) / 2) + 1
        #     elif i == 180:
        #         angles = -1
        #         means = mtr_means.shape[-1]
        #     plt.figure(figsize=(10, 5))
        #     plt.subplot(1, 2, 1)
        #     for idx in range(mtr_means.shape[0]):
        #         plt.scatter(bins[:angles], mtr_means[idx, :means], 2, label="FA thr = " + str(fa_thrs[idx]))
        #     plt.xlabel('Angle between e1 and B0 field (degrees)')
        #     plt.ylabel('MTR mean')
        #     plt.title('MTR vs Angle')
        #     plt.subplot(1, 2, 2)
        #     for idx in range(ihmtr_means.shape[0]):
        #         plt.scatter(bins[:angles], ihmtr_means[idx, :means], 2, label="FA thr = " + str(fa_thrs[idx]))
        #     plt.xlabel('Angle between e1 and B0 field (degrees)')
        #     plt.ylabel('ihMTR mean')
        #     plt.title('ihMTR vs Angle')
        #     plt.legend()
        #     plt.tight_layout()
        #     # plt.show()
        #     fig_name = "plot_all_FA_thr_" + str(j) + "_degrees_bins_NuFo_" + str(l) + "_" + str(i) + "_degrees_range" + ".png"
        #     fig_path = sub_ses_dir / fig_name
        #     plt.savefig(fig_path)
        #     plt.close()