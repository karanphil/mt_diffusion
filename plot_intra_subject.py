import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# This code actually plots all the results given as input in one graph. 
# It can be used for intrasubject or other things such as FA thr comparison.
# It should be run in the output folder.

argv = sys.argv # First input should be the output file and the rest should be the txt files.

def plot_init():
    plt.rcParams["font.family"] = "serif"
    plt.rcParams['font.serif'] = 'Helvetica'
    plt.style.use('seaborn-notebook')
    plt.rcParams['axes.grid'] = False
    plt.rcParams['grid.color'] = "darkgrey"
    plt.rcParams['grid.linewidth'] = 1
    plt.rcParams['grid.linestyle'] = "-"
    plt.rcParams['grid.alpha'] = "0.5"
    plt.rcParams['figure.figsize'] = (10.0, 5.0)
    plt.rcParams['font.size'] = 25
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.linewidth'] =2
    plt.rcParams['lines.linewidth']=2
    plt.rcParams['lines.markersize']=10
    plt.rcParams['text.latex.unicode']=True
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{amssymb}', r"\usepackage{amstext}"]
    # plt.rcParams['mathtext.default']='regular'

output = argv[1]

txt_files = []
for arg in argv[2:]:
    txt_files.append(arg)

temp_dims = np.loadtxt(txt_files[0], skiprows=1).shape
results = np.zeros(((len(txt_files),) + temp_dims))
for i, file in enumerate(txt_files):
    results[i] = np.loadtxt(file, skiprows=1)

bins = np.zeros((results[:, :, 0].shape[1] + 1))
bins[:180] = results[0, :, 0]
bins[-1] = results[0, -1, 1]

mtr_means = results[:, :, 2]
ihmtr_means = results[:, :, 3]

for i in [180]: # range of the angle bins to visualize
    # Plot the results
    if i == 90:
        angles = int((len(bins) - 1) / 2) + 1
        means = int((len(bins) - 1) / 2) + 1
    elif i == 180:
        angles = -1
        means = mtr_means.shape[-1]
    plot_init()
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    for idx in range(mtr_means.shape[0]):
        plt.scatter(bins[:angles], mtr_means[idx, :means], 2)
    plt.xlabel('Angle between e1 and B0 field (degrees)')
    plt.ylabel('MTR mean')
    plt.title('MTR vs Angle')
    plt.subplot(1, 2, 2)
    for idx in range(ihmtr_means.shape[0]):
        plt.scatter(bins[:angles], ihmtr_means[idx, :means], 2)
    plt.xlabel('Angle between e1 and B0 field (degrees)')
    plt.ylabel('ihMTR mean')
    plt.title('ihMTR vs Angle')
    # plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig(output) # does not take i into account!
    plt.close()

