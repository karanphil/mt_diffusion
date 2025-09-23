import argparse
from cmcrameri import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

from modules.io import plot_init


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('in_mtr',
                   help='Input MTR image.')
    
    p.add_argument('in_fixel_mtr',
                   help='Input fixel MTR image.')
    
    p.add_argument('out_image',
                   help='Path of the output image.')
    
    p.add_argument('reference',
                   help='Reference image.')
    
    p.add_argument('--mask', default=[],
                   help='Mask image.')
    
    p.add_argument('--thr', type=float, default=0.75)
    
    p.add_argument('--axis', default='axial', 
                   choices=['axial', 'sagittal', 'coronal'],
                   help='Axis to display the slice.')
    
    p.add_argument('--slice', type=int,
                   help='Index for where to slice the images.')

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    ref = nib.load(args.reference).get_fdata()

    data_shape = ref.shape

    mtr_img = nib.load(args.in_mtr)
    mtr = mtr_img.get_fdata().astype(np.float32)

    fixel_mtr_img = nib.load(args.in_fixel_mtr)
    fixel_mtr = fixel_mtr_img.get_fdata().astype(np.float32)

    if args.mask:
        mask = nib.load(args.mask).get_fdata()
        mask = (mask >= args.thr)
    else:
        mask = np.ones((data_shape))

    mtr *= mask
    mtr = np.ma.masked_where(mtr == 0, mtr)
    fixel_mtr *= mask
    fixel_mtr = np.ma.masked_where(fixel_mtr == 0, fixel_mtr)

    plot_init(dims=(20, 15), font_size=20)

    COLOR = 'white'
    mpl.rcParams['text.color'] = COLOR
    mpl.rcParams['axes.labelcolor'] = COLOR
    mpl.rcParams['xtick.color'] = COLOR
    mpl.rcParams['ytick.color'] = COLOR

    fig, ax = plt.subplots(1, 2, layout='constrained')

    if args.axis == 'axial':
        axis = 1
        if args.slice is not None:
            z_index = args.slice
        else:
            # find slice with highest sum of intensities
            z_index = np.argmax(np.sum(mtr, axis=(0, 1)))
        # x_index and y_index span across all shape
        x_index= slice(0, data_shape[0])
        y_index = slice(0, data_shape[1])
    elif args.axis == 'sagittal':
        axis = 1
        if args.slice is not None:
            x_index = args.slice
        else:
            # find slice with highest sum of intensities
            x_index = np.argmax(np.sum(mtr, axis=(1, 2)))
        y_index = slice(0, data_shape[1])
        z_index = slice(0, data_shape[2])
    else:  # coronal
        axis = 1
        if args.slice is not None:
            y_index = args.slice
        else:
            # find slice with highest sum of intensities
            y_index = np.argmax(np.sum(mtr, axis=(0, 2)))
        x_index= slice(0, data_shape[0])
        z_index = slice(0, data_shape[2])

    print(x_index, y_index, z_index)

    mtr_slice = np.flip(np.rot90(mtr[x_index, y_index, z_index]), axis=axis)
    fixel_mtr_slice = np.flip(np.rot90(fixel_mtr[x_index, y_index, z_index]),
                              axis=axis)

    colorbar = ax[0].imshow(mtr_slice, cmap=cm.navia, vmin=0, vmax=0.5,
                            interpolation='none')
    colorbar = ax[1].imshow(fixel_mtr_slice, cmap=cm.navia, vmin=0, vmax=0.5,
                            interpolation='none')

    # cb = fig.colorbar(colorbar, ax=ax[1], location='right', aspect=20, pad=0.1, fraction=0.01)
    # cb.outline.set_color('white')

    for i in range(ax.shape[0]):
        ax[i].set_axis_off()
        ax[i].autoscale(False)

    fig.get_layout_engine().set(h_pad=0.1, hspace=0.1) #, w_pad=0, wspace=0)
    # fig.tight_layout()
    # plt.show()
    plt.savefig(args.out_image, dpi=300, transparent=False, facecolor='black')


if __name__ == "__main__":
    main()
