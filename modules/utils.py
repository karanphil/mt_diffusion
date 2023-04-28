import numpy as np


def compute_single_fiber_averages(peaks, fa, wm_mask, affine, mtr=None,
                                  ihmtr=None, mtsat=None, ihmtsat=None,
                                  nufo=None, bin_width=1, fa_thr=0.5):
    peaks, _ = normalize_peaks(np.copy(peaks))
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


def compute_crossing_fibers_averages(peaks, wm_mask, affine, nufo, mtr=None,
                                     ihmtr=None, mtsat=None, ihmtsat=None,
                                     bin_width=10, norm_thr=0.7):
    peaks, peaks_norm = normalize_peaks(np.copy(peaks))

    # Find the direction of the B0 field
    rot = affine[0:3, 0:3]
    z_axis = np.array([0, 0, 1])
    b0_field = np.dot(rot.T, z_axis)

    bins = np.arange(0, 90 + bin_width, bin_width)
    # mid_bins = (bins[:-1] + bins[1:]) / 2.

    # Calculate the angle between e1 and B0 field
    cos_theta_f1 = np.dot(peaks[..., 0:3], b0_field)
    theta_f1 = np.arccos(cos_theta_f1) * 180 / np.pi
    cos_theta_f2 = np.dot(peaks[..., 3:6], b0_field)
    theta_f2 = np.arccos(cos_theta_f2) * 180 / np.pi

    mtr_means = np.zeros((len(bins) - 1, len(bins) - 1))
    ihmtr_means = np.zeros((len(bins) - 1, len(bins) - 1))
    mtsat_means = np.zeros((len(bins) - 1, len(bins) - 1))
    ihmtsat_means = np.zeros((len(bins) - 1, len(bins) - 1))
    nb_voxels = np.zeros((len(bins) - 1, len(bins) - 1))

    # Apply the WM mask
    wm_mask_bool = (wm_mask > 0.9) & (nufo == 2)
    norm_mask_bool = (peaks_norm[..., 0] > norm_thr) & (peaks_norm[..., 1] > norm_thr)

    for i in range(len(bins) - 1):
        angle_mask_0_90 = (theta_f1 >= bins[i]) & (theta_f1 < bins[i+1])
        angle_mask_90_180 = (180 - theta_f1 >= bins[i]) & (180 - theta_f1 < bins[i+1])
        angle_mask = angle_mask_0_90 | angle_mask_90_180
        mask_f1 = angle_mask
        for j in range(len(bins) - 1):
            angle_mask_0_90 = (theta_f2 >= bins[j]) & (theta_f2 < bins[j+1]) 
            angle_mask_90_180 = (180 - theta_f2 >= bins[j]) & (180 - theta_f2 < bins[j+1])
            angle_mask = angle_mask_0_90 | angle_mask_90_180
            mask_f2 = angle_mask
            mask = mask_f1 & mask_f2 & wm_mask_bool & norm_mask_bool
            nb_voxels[i, j] = np.sum(mask)
            if np.sum(mask) < 5:
                mtr_means[i, j] = None
                ihmtr_means[i, j] = None
                mtsat_means[i, j] = None
                ihmtsat_means[i, j] = None
            else:
                if mtr is not None:
                    mtr_means[i, j] = np.mean(mtr[mask])
                if ihmtr is not None:
                    ihmtr_means[i, j] = np.mean(ihmtr[mask])
                if mtsat is not None:
                    mtsat_means[i, j] = np.mean(mtsat[mask])
                if ihmtsat is not None:
                    ihmtsat_means[i, j] = np.mean(ihmtsat[mask])
    
    for i in range(len(bins) - 1):
        for j in range(i):
            mtr_means[i, j] = (mtr_means[i, j] + mtr_means[j, i]) / 2
            mtr_means[j, i] = mtr_means[i, j]
            ihmtr_means[i, j] = (ihmtr_means[i, j] + ihmtr_means[j, i]) / 2
            ihmtr_means[j, i] = ihmtr_means[i, j]
            mtsat_means[i, j] = (mtsat_means[i, j] + mtsat_means[j, i]) / 2
            mtsat_means[j, i] = mtsat_means[i, j]
            ihmtsat_means[i, j] = (ihmtsat_means[i, j] + ihmtsat_means[j, i]) / 2
            ihmtsat_means[j, i] = ihmtsat_means[i, j]
            nb_voxels[i, j] = nb_voxels[i, j] + nb_voxels[j, i]
            nb_voxels[j, i] = nb_voxels[i, j]

    return bins, mtr_means, ihmtr_means, mtsat_means, ihmtsat_means, nb_voxels


def extend_measure(bins, measure):
    new_bins = np.concatenate((np.flip(-bins[1:6]), bins, 180 - np.flip(bins[-6:-1])))
    new_measure = np.concatenate((np.flip(measure[1:6]), measure, np.flip(measure[-6:-1])))
    return new_bins, new_measure


def fit_single_fiber_results(bins, means, poly_order=8):
    new_bins, new_means = extend_measure(bins, means)
    mid_bins = (new_bins[:-1] + new_bins[1:]) / 2.
    not_nan = np.isfinite(new_means)
    fit = np.polyfit(mid_bins[not_nan], new_means[not_nan], poly_order)
    polynome = np.poly1d(fit)
    return polynome


def normalize_peaks(peaks):
    # Attention, il faut passer peaks en copie!!!
    nb_peaks = int(peaks.shape[-1] / 3)
    peaks_norm = np.zeros((peaks.shape[:-1] + (nb_peaks,)))
    for i in range(nb_peaks):
        peaks_norm[..., i] = np.linalg.norm(peaks[..., i * 3 : (i + 1) * 3], axis=-1)
        for j in range(3):
            peaks[..., i * 3 + j] /= peaks_norm[..., i]
    return peaks, peaks_norm


def compute_corrections(polynome, angle, fraction):
    bins = np.arange(0, 90 + 1, 1)
    max_poly = np.max(polynome(bins))
    correction = fraction * (max_poly - polynome(angle))
    return correction


def correct_measure(peaks, peak_values, measure, affine, wm_mask, polynome,
                    peak_frac_thr=0):
    wm_mask_bool = (wm_mask > 0.9)
    # peaks, peaks_norm = normalize_peaks(np.copy(peaks))
    # peaks_sum = np.sum(peaks_norm, axis=-1)
    # peaks_sum = np.repeat(peaks_sum.reshape(peaks_sum.shape + (1,)),
    #                       peaks_norm.shape[-1], axis=-1)
    # peaks_fraction = peaks_norm / peaks_sum
    max_peak_value = np.max(peak_values[wm_mask_bool])
    peaks_fraction = peak_values / max_peak_value
    
    # Find the direction of the B0 field
    rot = affine[0:3, 0:3]
    z_axis = np.array([0, 0, 1])
    b0_field = np.dot(rot.T, z_axis)

    bins = np.arange(0, 90 + 1, 1)
    # mid_bins = (bins[:-1] + bins[1:]) / 2.

    peaks_angles = np.empty((peaks_fraction.shape))
    peaks_angles[:] = np.nan
    corrections = np.zeros((peaks_fraction.shape))
    # Calculate the angle between e1 and B0 field for each peak
    for i in range(peaks_angles.shape[-1]):
        mask = wm_mask_bool & (peaks_fraction[..., i] > peak_frac_thr)
        cos_theta = np.dot(peaks[mask, i*3:(i+1)*3], b0_field)
        theta = np.arccos(cos_theta) * 180 / np.pi
        peaks_angles[mask, i] = np.abs(theta//90 * 90 - theta%90) % 180

        corrections[mask, i] = compute_corrections(polynome,
                                                   peaks_angles[mask, i],
                                                   peaks_fraction[mask, i])
    
    total_corrections = np.sum(corrections, axis=-1)

    return measure + total_corrections