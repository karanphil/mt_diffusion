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
    plt.rcParams['lines.markersize']=4
    # plt.rcParams['text.latex.unicode']=True
    # plt.rcParams['text.latex.preamble'] = [r'\usepackage{amssymb}', r"\usepackage{amstext}"]
    # plt.rcParams['mathtext.default']='regular'

e1_nifti = Path(argv[1])
fa_nifti = Path(argv[2])
wm_mask_nifti = Path(argv[3])
nufo_nifti = Path(argv[4])
sp_nifti = Path(argv[5])
sn_nifti = Path(argv[6])
sdualpn_nifti = Path(argv[7])
sdualnp_nifti = Path(argv[8])
ref_nifti = Path(argv[9])                     

# Load the data
sp_img = nib.load(str(sp_nifti))
sn_img = nib.load(str(sn_nifti))
sdualpn_img = nib.load(str(sdualpn_nifti))
sdualnp_img = nib.load(str(sdualnp_nifti))
fa_img = nib.load(str(fa_nifti))
wm_mask_img = nib.load(str(wm_mask_nifti))
e1_img = nib.load(str(e1_nifti))
nufo_img = nib.load(str(nufo_nifti))
ref_img = nib.load(str(ref_nifti))

sp_data = sp_img.get_fdata()
sn_data = sn_img.get_fdata()
sdualpn_data = sdualpn_img.get_fdata()
sdualnp_data = sdualnp_img.get_fdata()
fa_data = fa_img.get_fdata()
wm_mask_data = wm_mask_img.get_fdata()
e1_data = e1_img.get_fdata()
nufo_data = nufo_img.get_fdata()
ref_data = ref_img.get_fdata()

sub_ses_dir = Path(fa_nifti.parent.name)
if not sub_ses_dir.is_dir():
    sub_ses_dir.mkdir()

# bins_width = [1, 3, 5, 10]
# fa_thrs = [0.7, 0.6, 0.5]
bins_width = [1]
fa_thrs = [0.5]

# Find the direction of the B0 field
rot = fa_img.affine[0:3, 0:3]
z_axis = np.array([0, 0, 1])
b0_field = np.dot(rot.T, z_axis)

for j in bins_width: # width of the angle bins

    # Define the bins
    # bins = np.arange(0, 180, 5)
    bins = np.arange(0, 90 + j, j)
    mid_bins = (bins[:-1] + bins[1:]) / 2.

    # Calculate the angle between e1 and B0 field
    cos_theta = np.dot(e1_data[..., :3], b0_field)
    theta = np.arccos(cos_theta) * 180 / np.pi

    for l in [True]: # use NuFo or not
        sp_means = np.zeros((len(fa_thrs), len(bins) - 1))
        sn_means = np.zeros((len(fa_thrs), len(bins) - 1))
        spn_means = np.zeros((len(fa_thrs), len(bins) - 1))
        sdualpn_means = np.zeros((len(fa_thrs), len(bins) - 1))
        sdualnp_means = np.zeros((len(fa_thrs), len(bins) - 1))
        sdualpnnp_means = np.zeros((len(fa_thrs), len(bins) - 1))
        ref_means = np.zeros((len(fa_thrs), len(bins) - 1))
        nb_voxels = np.zeros((len(fa_thrs), len(bins) - 1))

        for k, fa_thr in enumerate(fa_thrs): # FA threshold

            # Apply the WM mask and FA threshold
            if l:
                wm_mask_bool = (wm_mask_data > 0.9) & (fa_data > fa_thr) & (nufo_data == 1) # & (rd_data > 5e-4) test rd thr # & (fa_data < 0.7) test upper fa thr
            else:
                wm_mask_bool = (wm_mask_data > 0.9) & (fa_data > fa_thr)

            for i in range(len(bins) - 1):
                angle_mask_0_90 = (theta >= bins[i]) & (theta < bins[i+1]) 
                angle_mask_90_180 = (180 - theta >= bins[i]) & (180 - theta < bins[i+1])
                angle_mask = angle_mask_0_90 | angle_mask_90_180
                mask = wm_mask_bool & angle_mask
                nb_voxels[k, i] = np.sum(mask)
                if np.sum(mask) < 5:
                    sp_means[k, i] = None
                    sn_means[k, i] = None
                    spn_means[k, i] = None
                    sdualpn_means[k, i] = None
                    sdualnp_means[k, i] = None
                    sdualpnnp_means[k, i] = None
                    ref_means[k, i] = None
                else:
                    sp_means[k, i] = np.mean(sp_data[mask])
                    sn_means[k, i] = np.mean(sn_data[mask])
                    spn_means[k, i] = np.mean((sp_data[mask] + sn_data[mask]) / 2)
                    sdualpn_means[k, i] = np.mean(sdualpn_data[mask])
                    sdualnp_means[k, i] = np.mean(sdualnp_data[mask])
                    sdualpnnp_means[k, i] = np.mean((sdualpn_data[mask] + sdualnp_data[mask]) / 2)
                    ref_means[k, i] = np.mean(ref_data[mask])

            # # Save the results to a text file
            # results = np.column_stack((bins[:-1], bins[1:], mtr_means[k], ihmtr_means[k], nb_voxels[k]))
            # txt_name = "results_" + str(j) + "_degrees_bins_" + str(fa_thr) + "_FA_thr_NuFo_" + str(l) + ".txt"
            # txt_path = sub_ses_dir / txt_name
            # if ratios:
            #     np.savetxt(str(txt_path), results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTR_mean\tihMTR_mean\tNb_voxels')
            # elif sats:
            #     np.savetxt(str(txt_path), results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTsat_mean\tihMTsat_mean\tNb_voxels')

            # Plot the results
            max_count = np.max(nb_voxels[k, :])
            norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
            plot_init()
            fig, (ax1, ax2, cax) = plt.subplots(1, 3, gridspec_kw={"width_ratios":[1,1, 0.05]})
            # ax1.scatter(mid_bins, sp_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='green', linewidths=1)
            # ax1.scatter(mid_bins, sn_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='blue', linewidths=1)
            ax1.scatter(mid_bins, spn_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='blue', linewidths=1)
            ax1.scatter(mid_bins, ref_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='green', linewidths=1)
            ax1.set_xlabel('Angle between e1 and B0 field (degrees)')
            ax1.set_ylabel('S+- mean')
            ax1.set_title('S+- vs Angle')
            # ax2.scatter(mid_bins, sdualpn_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='green', linewidths=1)
            # ax2.scatter(mid_bins, sdualnp_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='blue', linewidths=1)
            # colorbar = ax2.scatter(mid_bins, sdualpnnp_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='red', linewidths=1)
            colorbar = ax2.scatter(mid_bins, (ref_means[k, :] - spn_means[k, :]) / ref_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='red', linewidths=1)
            ax2.set_xlabel('Angle between e1 and B0 field (degrees)')
            ax2.set_ylabel('Sdual +-/-+ mean')
            ax2.set_title('Sdual +-/-+ vs Angle')
            fig.colorbar(colorbar, cax=cax, label="Voxel count")
            fig.tight_layout()
            # plt.show()
            fig_name = "plot_" + str(j) + "_degrees_bins_" + str(fa_thr) + "_FA_thr_NuFo_" + str(l) + "_Spn.png"
            fig_path = sub_ses_dir / fig_name
            plt.savefig(fig_path, dpi=300)
            plt.close()

            # Plot the results
            max_count = np.max(nb_voxels[k, :])
            norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
            plot_init()
            fig, (ax1, ax2, cax) = plt.subplots(1, 3, gridspec_kw={"width_ratios":[1,1, 0.05]})
            # ax1.scatter(mid_bins, sp_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='green', linewidths=1)
            # ax1.scatter(mid_bins, sn_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='blue', linewidths=1)
            ax1.scatter(mid_bins, sdualpnnp_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='blue', linewidths=1)
            ax1.scatter(mid_bins, ref_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='green', linewidths=1)
            ax1.set_xlabel('Angle between e1 and B0 field (degrees)')
            ax1.set_ylabel('S+- mean')
            ax1.set_title('S+- vs Angle')
            # ax2.scatter(mid_bins, sdualpn_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='green', linewidths=1)
            # ax2.scatter(mid_bins, sdualnp_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='blue', linewidths=1)
            # colorbar = ax2.scatter(mid_bins, sdualpnnp_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='red', linewidths=1)
            colorbar = ax2.scatter(mid_bins, (spn_means[k, :] - sdualpnnp_means[k, :]) / ref_means[k, :], c=nb_voxels[k, :], cmap='Greys', norm=norm, edgecolors='red', linewidths=1)
            ax2.set_xlabel('Angle between e1 and B0 field (degrees)')
            ax2.set_ylabel('Sdual +-/-+ mean')
            ax2.set_title('Sdual +-/-+ vs Angle')
            fig.colorbar(colorbar, cax=cax, label="Voxel count")
            fig.tight_layout()
            # plt.show()
            fig_name = "plot_" + str(j) + "_degrees_bins_" + str(fa_thr) + "_FA_thr_NuFo_" + str(l) + "_Sdualpn.png"
            fig_path = sub_ses_dir / fig_name
            plt.savefig(fig_path, dpi=300)
            plt.close()