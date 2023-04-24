import numpy as np


def compute_single_fiber_averages(peaks, fa, wm_mask, affine, mtr=None,
                                  ihmtr=None, mtsat=None, ihmtsat=None,
                                  nufo=None, bin_width=10, fa_thr=0.5):
    # Find the direction of the B0 field
    rot = affine[0:3, 0:3]
    z_axis = np.array([0, 0, 1])
    b0_field = np.dot(rot.T, z_axis)

    bins = np.arange(0, 90 + bin_width, bin_width)
    # mid_bins = (bins[:-1] + bins[1:]) / 2.

    # Calculate the angle between e1 and B0 field
    cos_theta = np.dot(peaks[..., :3], b0_field)
    theta = np.arccos(cos_theta) * 180 / np.pi

    mtr_means = np.zeros((len(bins) - 1))
    ihmtr_means = np.zeros((len(bins) - 1))
    mtsat_means = np.zeros((len(bins) - 1))
    ihmtsat_means = np.zeros((len(bins) - 1))
    nb_voxels = np.zeros((len(bins) - 1))

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
        nb_voxels[i] = np.sum(mask)
        if np.sum(mask) < 5:
            mtr_means[i] = None
            ihmtr_means[i] = None
            mtsat_means[i] = None
            ihmtsat_means[i] = None
        else:
            if mtr is not None:
                mtr_means[i] = np.mean(mtr[mask])
            if ihmtr is not None:
                ihmtr_means[i] = np.mean(ihmtr[mask])
            if mtsat is not None:
                mtsat_means[i] = np.mean(mtsat[mask])
            if ihmtsat is not None:
                ihmtsat_means[i] = np.mean(ihmtsat[mask])

    return bins, mtr_means, ihmtr_means, mtsat_means, ihmtsat_means, nb_voxels
