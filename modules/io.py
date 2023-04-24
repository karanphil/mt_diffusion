import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def save_txt(bins, mt_means, ihmt_means, nb_voxels, output_name, input_dtype="ratios"):
    # Save the results to a text file
    results = np.column_stack((bins[:-1], bins[1:], mt_means, ihmt_means, nb_voxels))
    if input_dtype=="ratios":
        np.savetxt(output_name, results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTR_mean\tihMTR_mean\tNb_voxels')
    elif input_dtype=="sats":
        np.savetxt(output_name, results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTsat_mean\tihMTsat_mean\tNb_voxels')


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


def plot_means(bins, mt_means, ihmt_means, nb_voxels, output_name, input_dtype="ratios"):
    max_count = np.max(nb_voxels)
    norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
    mid_bins = (bins[:-1] + bins[1:]) / 2.
    plot_init()
    fig, (ax1, ax2, cax) = plt.subplots(1, 3, gridspec_kw={"width_ratios":[1,1, 0.05]})
    ax1.scatter(mid_bins, mt_means, c=nb_voxels, cmap='Greys', norm=norm, edgecolors='black', linewidths=1)
    ax1.set_xlabel('Angle between e1 and B0 field (degrees)')
    if input_dtype=="ratios":
        ax1.set_ylabel('MTR mean')
        ax1.set_title('MTR vs Angle')
    elif input_dtype=="sats":
        ax1.set_ylabel('MTsat mean')
        ax1.set_title('MTsat vs Angle')
    colorbar = ax2.scatter(mid_bins, ihmt_means, c=nb_voxels, cmap='Greys', norm=norm, edgecolors='black', linewidths=1)
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