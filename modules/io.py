import nibabel as nib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path


def save_txt(bins, mt_means, ihmt_means, nb_voxels, output_name,
             input_dtype="ratios"):
    # Save the results to a text file
    results = np.column_stack((bins[:-1], bins[1:], mt_means, ihmt_means,
                               nb_voxels))
    if input_dtype=="ratios":
        np.savetxt(output_name, results, fmt='%10.5f', delimiter='\t',
                   header='Angle_min\tAngle_max\tMTR_mean\tihMTR_mean\tNb_voxels')
    elif input_dtype=="sats":
        np.savetxt(output_name, results, fmt='%10.5f', delimiter='\t',
                   header='Angle_min\tAngle_max\tMTsat_mean\tihMTsat_mean\tNb_voxels')


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


def plot_means(bins, mt_means, ihmt_means, nb_voxels, output_name,
               mt_cr_means=None, ihmt_cr_means=None,
               mt_poly=None, ihmt_poly=None, input_dtype="ratios"):
    max_count = np.max(nb_voxels)
    norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
    mid_bins = (bins[:-1] + bins[1:]) / 2.
    highres_bins = np.arange(0, 90 + 1, 0.5)
    plot_init()
    fig, (ax1, ax2, cax) = plt.subplots(1, 3,
                                        gridspec_kw={"width_ratios":[1,1, 0.05]})
    ax1.scatter(mid_bins, mt_means, c=nb_voxels, cmap='Greys', norm=norm,
                edgecolors='black', linewidths=1)
    if mt_cr_means is not None:
        ax1.scatter(mid_bins, mt_cr_means, c=nb_voxels, cmap='Greys',
                    norm=norm, edgecolors='blue', linewidths=1)
    if mt_poly:
        ax1.plot(highres_bins, mt_poly(highres_bins), "--g")
    ax1.set_xlabel('Angle between e1 and B0 field (degrees)')
    if input_dtype=="ratios":
        ax1.set_ylabel('MTR mean')
        ax1.set_title('MTR vs Angle')
    elif input_dtype=="sats":
        ax1.set_ylabel('MTsat mean')
        ax1.set_title('MTsat vs Angle')
    colorbar = ax2.scatter(mid_bins, ihmt_means, c=nb_voxels, cmap='Greys',
                           norm=norm, edgecolors='black', linewidths=1)
    if ihmt_cr_means is not None:
        ax2.scatter(mid_bins, ihmt_cr_means, c=nb_voxels, cmap='Greys',
                    norm=norm, edgecolors='blue', linewidths=1)
    if ihmt_poly:
        ax2.plot(highres_bins, ihmt_poly(highres_bins), "--g")
    ax2.set_xlabel('Angle between e1 and B0 field (degrees)')
    if input_dtype=="ratios":
        ax2.set_ylabel('ihMTR mean')
        ax2.set_title('ihMTR vs Angle')
    elif input_dtype=="sats":
        ax2.set_ylabel('ihMTsat mean')
        ax2.set_title('ihMTsat vs Angle')
    fig.colorbar(colorbar, cax=cax, label="Voxel count")
    fig.tight_layout()
    # plt.show()
    plt.savefig(output_name, dpi=300)
    plt.close()


def save_masks_by_angle_bins(peaks, fa, wm_mask, affine, output_path,
                             nufo=None, bin_width=10, fa_thr=0.5):
    # Find the direction of the B0 field
    rot = affine[0:3, 0:3]
    z_axis = np.array([0, 0, 1])
    b0_field = np.dot(rot.T, z_axis)

    # Define the bins
    bins = np.arange(0, 90 + bin_width, bin_width)

    # Calculate the angle between e1 and B0 field
    cos_theta = np.dot(peaks[..., :3], b0_field)
    theta = np.arccos(cos_theta) * 180 / np.pi

    # Apply the WM mask and FA threshold
    if nufo is not None:
        wm_mask_bool = (wm_mask > 0.9) & (fa > fa_thr) & (nufo == 1)
    else:
        wm_mask_bool = (wm_mask > 0.9) & (fa > fa_thr)
    for i in range(len(bins) - 1):
        angle_mask_0_90 = (theta >= bins[i]) & (theta < bins[i+1]) 
        angle_mask_90_180 = (180 - theta >= bins[i]) & (180 - theta < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask = wm_mask_bool & angle_mask
        mask_name = "single_fiber_mask_" + str(bins[i]) + "_to_" + str(bins[i+1]) \
            + "_degrees_angle_bin_" + str(fa_thr) + "_fa_thr_" + str(bin_width) + "_bin_width.nii.gz"
        mask_path = output_path / "masks" / mask_name
        nib.save(nib.Nifti1Image(mask.astype(np.uint8), affine), mask_path)


def plot_measure_mean(bins, peak_values, output_name):
    mid_bins = (bins[:-1] + bins[1:]) / 2.
    plot_init()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(mid_bins, peak_values, linewidths=1)
    ax.set_xlabel('Angle between e1 and B0 field (degrees)')
    ax.set_ylabel('Measure mean')
    ax.set_title('Measure mean vs Angle')
    fig.tight_layout()
    # plt.show()
    plt.savefig(output_name, dpi=300)
    plt.close()