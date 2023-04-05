import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import sys

# This code actually plots all the results given as input in one graph. 
# It can be used for intrasubject or other things such as FA thr comparison.
# It should be run in the output intrasujbect folder.

argv = sys.argv # First input should be the output filename and the rest should be the txt files.

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

txt_files = []
for arg in argv[2:7]:
    txt_files.append(arg)

temp_dims = np.loadtxt(txt_files[0], skiprows=1).shape
results = np.zeros(((len(txt_files),) + temp_dims))
for i, file in enumerate(txt_files):
    results[i] = np.loadtxt(file, skiprows=1)

bins = np.zeros((results[:, :, 0].shape[1] + 1))
bins[:90] = results[0, :, 0]
bins[-1] = results[0, -1, 1]

mtr_means = results[:, :, 2]
ihmtr_means = results[:, :, 3]
nb_voxels = results[:, :, 4]

# Plot the results
max_count = np.max(nb_voxels[:, :])
norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
plot_init()
fig, (ax1, ax2, cax) = plt.subplots(1, 3, gridspec_kw={"width_ratios":[1,1, 0.05]})
for idx in range(mtr_means.shape[0]):
    ax1.scatter(bins[:-1], mtr_means[idx, :], c=nb_voxels[idx, :], cmap='Greys', norm=norm, edgecolors="C" + str(idx), linewidths=1, label="Session " + str(idx + 1))
ax1.set_xlabel('Angle between e1 and B0 field (degrees)')
ax1.set_ylabel('MTR mean')
ax1.set_title('MTR vs Angle')
for idx in range(ihmtr_means.shape[0]):
    colorbar = ax2.scatter(bins[:-1], ihmtr_means[idx, :], c=nb_voxels[idx, :], cmap='Greys', norm=norm, edgecolors="C" + str(idx), linewidths=1, label="Session " + str(idx + 1))
ax2.set_xlabel('Angle between e1 and B0 field (degrees)')
ax2.set_ylabel('ihMTR mean')
ax2.set_title('ihMTR vs Angle')
ax2.legend()
fig.colorbar(colorbar, cax=cax, label="Voxel count")
fig.tight_layout()
# plt.show()
output = output_name + "_" + str(i) + "_degrees_range" + ".png"
plt.savefig(output, dpi=300)
plt.close()