import nibabel as nib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path


def compute_peaks_fraction(peak_values):
    peak_values_sum = np.sum(peak_values, axis=-1)
    peak_values_sum = np.repeat(peak_values_sum.reshape(peak_values_sum.shape + (1,)),
                                peak_values.shape[-1], axis=-1)
    peaks_fraction = peak_values / peak_values_sum
    return peaks_fraction


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
    plt.rcParams['figure.figsize'] = (12.0, 5.0)
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
               mt_poly=None, ihmt_poly=None, input_dtype="ratios",
               labels=[None, None, None]):
    max_count = np.max(nb_voxels)
    norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
    mid_bins = (bins[:-1] + bins[1:]) / 2.
    highres_bins = np.arange(0, 90 + 1, 0.5)
    plot_init()
    fig, (ax1, ax2, cax) = plt.subplots(1, 3,
                                        gridspec_kw={"width_ratios":[1,1, 0.05]})
    ax1.scatter(mid_bins, mt_means, c=nb_voxels, cmap='Greys', norm=norm,
                edgecolors="C0", linewidths=1)
    if mt_cr_means is not None:
        ax1.scatter(mid_bins, mt_cr_means, c=nb_voxels, cmap='Greys',
                    norm=norm, edgecolors="C0", linewidths=1, marker="s")
    if mt_poly:
        ax1.plot(highres_bins, mt_poly(highres_bins), "--", color="C0")
    ax1.set_xlabel(r'$\theta_a$')
    ax1.set_xlim(0, 90)
    if input_dtype=="ratios":
        ax1.set_ylabel('MTR mean')
        # ax1.set_title('MTR vs Angle')
    elif input_dtype=="sats":
        ax1.set_ylabel('MTsat mean')
        # ax1.set_title('MTsat vs Angle')
    colorbar = ax2.scatter(mid_bins, ihmt_means, c=nb_voxels, cmap='Greys',
                           norm=norm, edgecolors="C0", linewidths=1, label=labels[0])
    if ihmt_cr_means is not None:
        ax2.scatter(mid_bins, ihmt_cr_means, c=nb_voxels, cmap='Greys',
                    norm=norm, edgecolors="C0", linewidths=1, marker="s", label=labels[1])
    if ihmt_poly:
        ax2.plot(highres_bins, ihmt_poly(highres_bins), "--", color="C0", label=labels[2])
    ax2.set_xlabel(r'$\theta_a$')
    ax2.set_xlim(0, 90)
    if input_dtype=="ratios":
        ax2.set_ylabel('ihMTR mean')
        # ax2.set_title('ihMTR vs Angle')
    elif input_dtype=="sats":
        ax2.set_ylabel('ihMTsat mean')
        # ax2.set_title('ihMTsat vs Angle')
    if labels[0] is not None:
        ax2.legend(loc=5)
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
        wm_mask_bool = (wm_mask >= 0.9) & (fa > fa_thr) & (nufo == 1)
    else:
        wm_mask_bool = (wm_mask >= 0.9) & (fa > fa_thr)
    for i in range(len(bins) - 1):
        angle_mask_0_90 = (theta >= bins[i]) & (theta < bins[i+1]) 
        angle_mask_90_180 = (180 - theta >= bins[i]) & (180 - theta < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask = wm_mask_bool & angle_mask
        mask_name = "single_fiber_mask_" + str(bins[i]) + "_to_" + str(bins[i+1]) \
            + "_degrees_angle_bin_" + str(fa_thr) + "_fa_thr_" + str(bin_width) + "_bin_width.nii.gz"
        mask_path = output_path / "masks" / mask_name
        nib.save(nib.Nifti1Image(mask.astype(np.uint8), affine), mask_path)


def save_angle_map(peaks, fa, wm_mask, affine, output_path, fodf_peaks, peak_values,
                    nufo, bin_width=1, fa_thr=0.5):
    # Find the direction of the B0 field
    rot = affine[0:3, 0:3]
    z_axis = np.array([0, 0, 1])
    b0_field = np.dot(rot.T, z_axis)

    # Define the bins
    bins = np.arange(0, 90 + bin_width, bin_width)

    # Calculate the angle between e1 and B0 field
    cos_theta = np.dot(peaks[..., :3], b0_field)
    theta = np.arccos(cos_theta) * 180 / np.pi

    peaks_fraction = compute_peaks_fraction(peak_values)

    cos_theta_f1 = np.dot(fodf_peaks[..., 0:3], b0_field)
    theta_f1 = np.arccos(cos_theta_f1) * 180 / np.pi
    cos_theta_f2 = np.dot(fodf_peaks[..., 3:6], b0_field)
    theta_f2 = np.arccos(cos_theta_f2) * 180 / np.pi
    cos_theta_f3 = np.dot(fodf_peaks[..., 6:9], b0_field)
    theta_f3 = np.arccos(cos_theta_f3) * 180 / np.pi

    peak_1 = np.zeros(wm_mask.shape)
    peak_2 = np.zeros(wm_mask.shape)
    peak_3 = np.zeros(wm_mask.shape)

    # Apply the WM mask and FA threshold
    wm_mask_bool = (wm_mask > 0.9) & (fa > fa_thr) & (nufo == 1)
    for i in range(len(bins) - 1):
        angle_mask_0_90 = (theta >= bins[i]) & (theta < bins[i+1]) 
        angle_mask_90_180 = (180 - theta >= bins[i]) & (180 - theta < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask = wm_mask_bool & angle_mask
        peak_1[mask] = (bins[i] + bins[i + 1]) /2.
    
    peak_1_sf = np.copy(peak_1)

    wm_mask_bool = (wm_mask > 0.9) & (nufo == 2)
    fraction_mask_bool = (peaks_fraction[..., 0] >= 0.5) & (peaks_fraction[..., 0] < 0.9)
    for i in range(len(bins) - 1):
        angle_mask_0_90 = (theta_f1 >= bins[i]) & (theta_f1 < bins[i+1])
        angle_mask_90_180 = (180 - theta_f1 >= bins[i]) & (180 - theta_f1 < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask_f1 = angle_mask & fraction_mask_bool & wm_mask_bool
        peak_1[mask_f1] = (bins[i] + bins[i + 1]) /2.

        angle_mask_0_90 = (theta_f2 >= bins[i]) & (theta_f2 < bins[i+1]) 
        angle_mask_90_180 = (180 - theta_f2 >= bins[i]) & (180 - theta_f2 < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask_f2 = angle_mask & fraction_mask_bool & wm_mask_bool
        peak_2[mask_f2] = (bins[i] + bins[i + 1]) /2.

    wm_mask_bool = (wm_mask > 0.9) & (nufo == 3)
    fraction_mask_bool = (peaks_fraction[..., 0] >= 0.33) & (peaks_fraction[..., 0] < 0.8)
    for i in range(len(bins) - 1):
        angle_mask_0_90 = (theta_f1 >= bins[i]) & (theta_f1 < bins[i+1])
        angle_mask_90_180 = (180 - theta_f1 >= bins[i]) & (180 - theta_f1 < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask_f1 = angle_mask & wm_mask_bool & fraction_mask_bool
        peak_1[mask_f1] = (bins[i] + bins[i + 1]) /2.

        angle_mask_0_90 = (theta_f2 >= bins[i]) & (theta_f2 < bins[i+1]) 
        angle_mask_90_180 = (180 - theta_f2 >= bins[i]) & (180 - theta_f2 < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask_f2 = angle_mask & wm_mask_bool & fraction_mask_bool
        peak_2[mask_f2] = (bins[i] + bins[i + 1]) /2.

        angle_mask_0_90 = (theta_f3 >= bins[i]) & (theta_f3 < bins[i+1]) 
        angle_mask_90_180 = (180 - theta_f3 >= bins[i]) & (180 - theta_f3 < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask_f3 = angle_mask & wm_mask_bool & fraction_mask_bool
        peak_3[mask_f3] = (bins[i] + bins[i + 1]) /2.

    map_1_name = "peak_1_sf_angles_map.nii.gz"
    map_1_path = output_path / "masks" / map_1_name
    nib.save(nib.Nifti1Image(peak_1_sf, affine), map_1_path)

    map_1_name = "peak_1_angles_map.nii.gz"
    map_1_path = output_path / "masks" / map_1_name
    nib.save(nib.Nifti1Image(peak_1, affine), map_1_path)

    map_2_name = "peak_2_angles_map.nii.gz"
    map_2_path = output_path / "masks" / map_2_name
    nib.save(nib.Nifti1Image(peak_2, affine), map_2_path)

    map_3_name = "peak_3_angles_map.nii.gz"
    map_3_path = output_path / "masks" / map_3_name
    nib.save(nib.Nifti1Image(peak_3, affine), map_3_path)


def plot_measure_mean(bins, peak_values, output_name):
    mid_bins = (bins[:-1] + bins[1:]) / 2.
    plot_init()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(mid_bins, peak_values, linewidths=1)
    ax.set_xlabel(r'$\theta_a$')
    ax.set_ylabel('Measure mean')
    # ax.set_title('Measure mean vs Angle')
    fig.tight_layout()
    # plt.show()
    plt.savefig(output_name, dpi=300)
    plt.close()


def plot_multiple_means(bins, mt_means, ihmt_means, nb_voxels, output_name,
                        labels=None, mt_cr_means=None, ihmt_cr_means=None,
                        input_dtype="ratios", legend_title=None,
                        mt_poly=None, ihmt_poly=None, delta_plot=False,
                        xlim=[0, 1.03], p_frac=None, mt_max=None,
                        mt_max_poly=None, leg_loc=1, markers="o"):
    max_count = np.max(nb_voxels)
    norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
    mid_bins = (bins[:-1] + bins[1:]) / 2.
    highres_bins = np.arange(0, 90 + 1, 0.5)
    plot_init()
    fig, (ax1, ax2, cax) = plt.subplots(1, 3,
                                        gridspec_kw={"width_ratios":[1,1, 0.05]})
    for i in range(mt_means.shape[0]):
        ax1.scatter(mid_bins, mt_means[i], c=nb_voxels[i], cmap='Greys', norm=norm,
                    linewidths=1, edgecolors="C" + str(i), marker=markers)
        if mt_cr_means is not None:
            ax1.scatter(mid_bins, mt_cr_means[i], c=nb_voxels[i], cmap='Greys',
                        norm=norm, edgecolors="C" + str(i), linewidths=1, marker="s")
        if mt_poly is not None:
            ax1.plot(highres_bins, mt_poly[i](highres_bins), "--", color="C" + str(i))
        if labels is not None:
            colorbar = ax2.scatter(mid_bins, ihmt_means[i], c=nb_voxels[i], cmap='Greys',
                                norm=norm, linewidths=1, label=labels[i], edgecolors="C" + str(i))
        else:
            colorbar = ax2.scatter(mid_bins, ihmt_means[i], c=nb_voxels[i], cmap='Greys',
                                norm=norm, linewidths=1, edgecolors="C" + str(i), marker=markers)
        if ihmt_cr_means is not None:
            ax2.scatter(mid_bins, ihmt_cr_means[i], c=nb_voxels[i], cmap='Greys',
                        norm=norm, edgecolors="C" + str(i), linewidths=1, marker="s")
        if ihmt_poly is not None:
            ax1.plot(highres_bins, ihmt_poly[i](highres_bins), "--", color="C" + str(i))
    ax1.set_xlabel(r'$\theta_a$')
    ax1.set_xlim(0, 90)
    if input_dtype=="ratios":
        ax1.set_ylabel('MTR mean')
        # ax1.set_title('MTR vs Angle')
    elif input_dtype=="sats":
        ax1.set_ylabel('MTsat mean')
        # ax1.set_title('MTsat vs Angle')

    if delta_plot:
        # this is an inset axes over the main axes
        highres_frac = np.arange(0, 1.01, 0.01)
        # ax = inset_axes(ax1,
        #                 bbox_to_anchor=[0.2, 0.2, 0.2, 0.2],
        #                 width="50%", # width = 40% of parent_bbox
        #                 height=1.0) # height : 1 inch
        #                 # loc=2)
        ax = plt.axes([0.13, 0.75, 0.16, 0.2])
        for i in range(len(p_frac) - 1):
            ax.scatter(p_frac[i], mt_max[i], color="C" + str(i), linewidths=1)
        ax.scatter(p_frac[-1], mt_max[-1], color="black", linewidths=1)
        ax.plot(highres_frac, mt_max_poly(highres_frac), "--", color="grey")
        ax.set_xlabel(r'Peak$_1$ fraction')
        ax.set_xlim(xlim[0], xlim[1])
        if input_dtype=="ratios":
            ax.set_ylabel(r'MTR $\delta m_{max}$')
        elif input_dtype=="sats":
            ax.set_ylabel(r'MTsat $\delta m_{max}$')
    
    ax2.set_xlabel(r'$\theta_a$')
    ax2.set_xlim(0, 90)
    if input_dtype=="ratios":
        ax2.set_ylabel('ihMTR mean')
        # ax2.set_title('ihMTR vs Angle')
    elif input_dtype=="sats":
        ax2.set_ylabel('ihMTsat mean')
        # ax2.set_title('ihMTsat vs Angle')
    if labels is not None:
        ax2.legend(loc=leg_loc)
    if legend_title is not None:
        ax2.get_legend().set_title(legend_title)
    fig.colorbar(colorbar, cax=cax, label="Voxel count")
    fig.tight_layout()
    # plt.show()
    plt.savefig(output_name, dpi=300)
    plt.close()


def plot_delta_m_max(bins, mt_means, mt_poly, output_name, ihmt_means=None,
                     ihmt_poly=None, input_dtype="ratios", xlim=[0, 1]):
    highres_bins = np.arange(0, 1.01, 0.01)
    plot_init()
    plt.rcParams['figure.figsize'] = (6.0, 5.0)
    if ihmt_means is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2,
                                            gridspec_kw={"width_ratios":[1,1]})
        ax1.scatter(bins, mt_means, color='black', linewidths=1)
        ax1.plot(highres_bins, mt_poly(highres_bins), "--g")
        ax1.set_xlabel(r'Peak$_1$ fraction')
        ax1.set_xlim(xlim[0], xlim[1])
        if input_dtype=="ratios":
            ax1.set_ylabel(r'MTR $\delta m_{max}$')
            # ax1.set_title('MTR vs Angle')
        elif input_dtype=="sats":
            ax1.set_ylabel(r'MTsat $\delta m_{max}$')
            # ax1.set_title('MTsat vs Angle')
        ax2.scatter(bins, ihmt_means, color='black', linewidths=1)
        ax2.plot(highres_bins, ihmt_poly(highres_bins), "--g")
        ax2.set_xlabel(r'Peak$_1$ fraction')
        ax2.set_xlim(xlim[0], xlim[1])
        if input_dtype=="ratios":
            ax2.set_ylabel(r'ihMTR $\delta m_{max}$')
            # ax2.set_title('ihMTR vs Angle')
        elif input_dtype=="sats":
            ax2.set_ylabel(r'ihMTsat $\delta m_{max}$')
            # ax2.set_title('ihMTsat vs Angle')
        fig.tight_layout()
        # plt.show()
        plt.savefig(output_name, dpi=300)
        plt.close()
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(bins, mt_means, color='black', linewidths=1)
        ax.plot(highres_bins, mt_poly(highres_bins), "--g")
        ax.set_xlabel(r'Peak$_1$ fraction')
        ax.set_xlim(xlim[0], xlim[1])
        if input_dtype=="ratios":
            ax.set_ylabel(r'MTR $\delta m_{max}$')
        elif input_dtype=="sats":
            ax.set_ylabel(r'MTsat $\delta m_{max}$')
        fig.tight_layout()
        # plt.show()
        plt.savefig(output_name, dpi=300)
        plt.close()


def plot_different_bins_means(bins, mt_means, ihmt_means, nb_voxels, output_name,
                              labels, input_dtype="ratios", legend_title=None,
                              mt_poly=None, ihmt_poly=None, mt_cr_means=None,
                              ihmt_cr_means=None, ax1_legend=False):
    max_count = 0
    for i in range(len(nb_voxels)):
        if np.max(nb_voxels[i]) > max_count:
            max_count = np.max(nb_voxels[i])
    norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
    highres_bins = np.arange(0, 90 + 1, 0.5)
    plot_init()
    fig, (ax1, ax2, cax) = plt.subplots(1, 3,
                                        gridspec_kw={"width_ratios":[1,1, 0.05]})
    mid_bins_1 = (bins[0][:-1] + bins[0][1:]) / 2.
    for i in range(len(mt_means)):
        mid_bins = (bins[i][:-1] + bins[i][1:]) / 2.
        ax1.scatter(mid_bins, mt_means[i], c=nb_voxels[i], cmap='Greys', norm=norm,
                    linewidths=1, label=labels[i], edgecolors="C" + str(i))
        if mt_cr_means is not None:
            ax1.scatter(mid_bins_1, mt_cr_means[i], c=nb_voxels[0], cmap='Greys',
                        norm=norm, edgecolors="C" + str(i), linewidths=1, marker="s")
        if mt_poly is not None:
            ax1.plot(highres_bins, mt_poly[i](highres_bins), "--", color="C" + str(i))
        colorbar = ax2.scatter(mid_bins, ihmt_means[i], c=nb_voxels[i], cmap='Greys',
                               norm=norm, linewidths=1, label=labels[i], edgecolors="C" + str(i))
        if ihmt_cr_means is not None:
            ax2.scatter(mid_bins_1, ihmt_cr_means[i], c=nb_voxels[0], cmap='Greys',
                        norm=norm, edgecolors="C" + str(i), linewidths=1, marker="s")
        if ihmt_poly is not None:
            ax2.plot(highres_bins, ihmt_poly[i](highres_bins), "--", color="C" + str(i))
    ax1.set_xlabel(r'$\theta_a$')
    ax1.set_xlim(0, 90)
    if input_dtype=="ratios":
        ax1.set_ylabel('MTR mean')
        # ax1.set_title('MTR vs Angle')
    elif input_dtype=="sats":
        ax1.set_ylabel('MTsat mean')
        # ax1.set_title('MTsat vs Angle')
    ax2.set_xlabel(r'$\theta_a$')
    ax2.set_xlim(0, 90)
    if input_dtype=="ratios":
        ax2.set_ylabel('ihMTR mean')
        # ax2.set_title('ihMTR vs Angle')
    elif input_dtype=="sats":
        ax2.set_ylabel('ihMTsat mean')
        # ax2.set_title('ihMTsat vs Angle')
    if ax1_legend:
        ax1.legend(loc=4)
        if legend_title is not None:
            ax1.get_legend().set_title(legend_title)
    else:
        ax2.legend()
        if legend_title is not None:
            ax2.get_legend().set_title(legend_title)
    fig.colorbar(colorbar, cax=cax, label="Voxel count")
    fig.tight_layout()
    # plt.show()
    plt.savefig(output_name, dpi=300)
    plt.close()


def plot_3d_means(bins, means, base_name, input_dtype="MTR"):
    mid_bins = (bins[:-1] + bins[1:]) / 2.
    plot_init()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(mid_bins, mid_bins)
    ax.plot_surface(X, Y, means, cmap="jet")
    ax.set_xlabel(r'$\theta_{a1}$')
    ax.set_ylabel(r'$\theta_{a2}$')
    ax.set_zlabel(input_dtype + ' mean')
    fig.tight_layout()
    # plt.show()
    if input_dtype=="MTR" or input_dtype=="MTsat":
        views = np.array([[30, -135], [30, -45], [10, -90], [10, 0]])
    elif input_dtype=="ihMTR" or input_dtype=="ihMTsat":
        views = np.array([[30, 45], [30, -45], [10, -90], [10, 0]])
    for v, view in enumerate(views[:]):
        output_name = base_name + "view_" + str(v) + ".png"
        ax.view_init(view[0], view[1])
        plt.savefig(output_name, dpi=300)
    plt.close()