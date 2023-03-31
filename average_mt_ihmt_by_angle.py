import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt


def main(argv):
    # Parse command line arguments
    if len(argv) != 6:
        print('Usage: python average_mt_ihmt_by_angle.py e1_nifti fa_nifti mtr_nifti ihmtr_nifti wm_mask_nifti')
        return
    e1_nifti = argv[1]
    fa_nifti = argv[2]
    mtr_nifti = argv[3]
    ihmtr_nifti = argv[4]
    wm_mask_nifti = argv[5]

    # Load the data
    mtr_img = nib.load(mtr_nifti)
    ihmtr_img = nib.load(ihmtr_nifti)
    fa_img = nib.load(fa_nifti)
    wm_mask_img = nib.load(wm_mask_nifti)
    e1_img = nib.load(e1_nifti)

    mtr_data = mtr_img.get_fdata()
    ihmtr_data = ihmtr_img.get_fdata()
    fa_data = fa_img.get_fdata()
    wm_mask_data = wm_mask_img.get_fdata()
    e1_data = e1_img.get_fdata()

    # Define the bins
    bins = np.arange(0, 180, 5)

    # Calculate the angle between e1 and B0 field
    b0_field = np.array([0, 0, 1])
    cos_theta = np.dot(e1_data[..., :3], b0_field)
    theta = np.arccos(cos_theta) * 180 / np.pi

    # Apply the WM mask and FA threshold
    wm_mask_bool = (wm_mask_data > 0) & (fa_data > 0.5)

    # Calculate the mean MTR and ihMTR for each bin
    mtr_means = np.zeros(len(bins) - 1)
    ihmtr_means = np.zeros(len(bins) - 1)

    for i in range(len(bins) - 1):
        angle_mask = (theta >= bins[i]) & (theta < bins[i+1])
        mask = wm_mask_bool & angle_mask
        mtr_means[i] = np.mean(mtr_data[mask])
        ihmtr_means[i] = np.mean(ihmtr_data[mask])

    # Save the results to a text file
    results = np.column_stack((bins[:-1], bins[1:], mtr_means, ihmtr_means))
    np.savetxt('angle_means.txt', results, fmt='%10.5f', delimiter='\t', header='Angle_min\tAngle_max\tMTR_mean\tihMTR_mean')

    # Plot the results
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(bins[:-1], mtr_means)
    plt.xlabel('Angle between e1 and B0 field (degrees)')
    plt.ylabel('MTR mean')
    plt.title('MTR vs Angle')
    plt.subplot(1, 2, 2)
    plt.plot(bins[:-1], ihmtr_means)
    plt.xlabel('Angle between e1 and B0 field (degrees)')
    plt.ylabel('ihMTR mean')
    plt.title('ihMTR vs Angle')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()