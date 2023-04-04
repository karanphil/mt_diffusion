import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from pathlib import Path
import sys

argv = sys.argv # First input should be the output file and the rest should be the txt files.

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

for i in [90, 180]: # range of the angle bins to visualize
    # Plot the results
    if i == 90:
        angles = int((len(bins) - 1) / 2) + 1
        means = int((len(bins) - 1) / 2) + 1
    elif i == 180:
        angles = -1
        means = mtr_means.shape[-1]
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
    plt.savefig(output)
    plt.close()

