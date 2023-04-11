from scipy.optimize import curve_fit
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys


def mt_curve(theta, a, b):
    theta_rad = np.deg2rad(theta)
    return a  * (3 * np.cos(theta_rad)**2 - 1) + b

# This code fits the MT and ihMT curves.
# It should be run in the output folder.

argv = sys.argv # First input should be the output file and the rest should be the subjects. Second and third should be output_choice (rations or sats) and bin_width.

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

output_name = argv[1]
output_choice = argv[2]
input_results = argv[3]

# results_name = "/results_" + str(bin_width) + "_degrees_bins_0.5_FA_thr_NuFo_False.txt"

# Select MTR and ihMTR or MTsat and ihMTsat
if output_choice == "ratios":
    ratios = True
    sats = False
elif output_choice == "sats":
    ratios = False
    sats = True

results = np.loadtxt(input_results, skiprows=1)

bins = np.zeros((results[:, 0].shape[0] + 1))
bins[:results[:, 0].shape[0]] = results[:, 0]
bins[-1] = results[-1, 1]
mid_bins = (bins[:-1] + bins[1:]) / 2.

mtr_means = results[:, 2]
ihmtr_means = results[:, 3]
nb_voxels = results[:, 4]

# mt_fit, mt_fit_cov = curve_fit(mt_curve, mid_bins, mtr_means)
mt_fit = np.polyfit(mid_bins, mtr_means, 10)
mt_p = np.poly1d(mt_fit)

ihmt_fit = np.polyfit(mid_bins, ihmtr_means, 10)
ihmt_p = np.poly1d(ihmt_fit)

# Plot the results
max_count = np.max(nb_voxels)
norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
plot_init()
fig, (ax1, ax2, cax) = plt.subplots(1, 3, gridspec_kw={"width_ratios":[1,1, 0.05]})
ax1.scatter(mid_bins, mtr_means, c=nb_voxels, cmap='Greys', norm=norm, edgecolors='black', linewidths=1)
# ax1.plot(mid_bins, mt_curve(mid_bins, 1, 0))
ax1.plot(mid_bins, mt_p(mid_bins), "--g")
ax1.set_xlabel('Angle between e1 and B0 field (degrees)')
if ratios:
    ax1.set_ylabel('MTR mean')
    ax1.set_title('MTR vs Angle')
elif sats:
    ax1.set_ylabel('MTsat mean')
    ax1.set_title('MTsat vs Angle')
colorbar = ax2.scatter(mid_bins, ihmtr_means, c=nb_voxels, cmap='Greys', norm=norm, edgecolors='black', linewidths=1)
ax2.plot(mid_bins, ihmt_p(mid_bins), "--g")
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
# fig_name = "plot_" + str(j) + "_degrees_bins_" + str(fa_thr) + "_FA_thr_NuFo_" + str(l) + ".png"
# fig_path = sub_ses_dir / fig_name
# plt.savefig(fig_path, dpi=300)
plt.close()