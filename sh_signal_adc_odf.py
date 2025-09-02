
import logging
import numpy as np
from dipy.sims.voxel import multi_tensor, multi_tensor_odf
from dipy.data import get_sphere
from dipy.core.sphere import disperse_charges, Sphere, HemiSphere
from dipy.core.gradients import gradient_table

from fury import window, actor

import nibabel as nib
from dipy.io.gradients import read_bvals_bvecs
from scilpy.gradients.bvec_bval_tools import (check_b0_threshold,
                                              is_normalized_bvecs,
                                              normalize_bvecs)

# Change this to change the size of your pop up windows.
WINDOW_SIZE = (600, 600)

def my_visu(sf, sphere, rot=True, norm=True, scale=True, title="Modeling"):
    ren = window.Scene(background=window.colors.white)
    
    sf_actor = actor.odf_slicer(sf[None, None, None, :], sphere=sphere, colormap='jet',
                                norm=norm)
    if rot :
        sf_actor.RotateX(90)

    ren.add(sf_actor)
    window.show(ren, title=title, size=WINDOW_SIZE)

    ren.rm(sf_actor)
    #window.rm_all(ren)

#
# generate gradient table 
#

img = nib.load("dwi_UNMASC_102855_zoom_in_septum.nii.gz")
data = img.get_fdata(dtype=np.float32)
affine = img.affine
bvals, bvecs = read_bvals_bvecs("bval.txt", "bvec.txt")

if not is_normalized_bvecs(bvecs):
    logging.warning('Your b-vectors do not seem normalized...')
    bvecs = normalize_bvecs(bvecs)

gtab = gradient_table(bvals, bvecs=bvecs)
gtab_vis = gradient_table(bvals[1:], bvecs=bvecs[1:,:])

print(bvecs)
print(bvals)

# on prend un voxel du signal dans le vert
signal = data[17, 25, 25, :]

signal_sph = np.zeros((signal.shape[0]-1)*2)
signal_sph[0:signal.shape[0]-1] = signal[1:]
signal_sph[signal.shape[0]-1:] = signal[1:]

print(signal_sph)
sphere = get_sphere('symmetric724')
sphere = sphere.subdivide(2)

print('Spherical Harmonics')
from dipy.reconst.shm import sh_to_sf_matrix
sh_order=8
vertices = gtab_vis.gradients
sph_gtab_vis = Sphere(xyz=np.vstack((vertices, -vertices)))
my_visu(signal_sph, sph_gtab_vis, False, False, False, title="Raw signal")

B, invB = sh_to_sf_matrix(sph_gtab_vis, sh_order)
B_highres, invB_highres = sh_to_sf_matrix(sphere, sh_order)

sh_signal = np.dot(invB.T, signal_sph)
signal_sphere = np.dot(sh_signal, B_highres)

my_visu(signal_sphere, sphere, False, False, False, title="Raw signal")


1/0





#
# Playing with the Apparent Diffusion Coefficient (ADC)
#
print('ADC')
adc = -1/bvalue*np.log(signal_sph)





# min-max normalization to enhance the spherical values
adc = (adc - adc.min()) / (adc.max()-adc.min())
if interactive :
    my_visu(adc, sph_gtab, True, True, True,
            title="ADC - min-max normalized")

    
#
# Playing with Spherical Harmonics
#

print('Building high resolution SH matrix of order 8')
sphere = get_sphere('symmetric724')
sphere = sphere.subdivide(2)
B_highres, invB_highres = sh_to_sf_matrix(sphere, 8)
print(B_highres.shape)

print('Estimate SH coefficients of the signal')
# S = B*sh_signal  => sh_signal = invB * S

print(sh_signal.shape)
print('Project back SH coefficient to the high resolution sphere')
signal_sphere = np.dot(sh_signal, B_highres)
print(signal_sphere.shape)

if interactive :
    my_visu(signal_sph, sph_gtab, False,
            title="Basic signal on the sphere")
    my_visu(signal_sphere, sphere, False,
            title="SH backprojected on High-res sphere")

print('SH coefficients of order 8')
print(sh_signal)
print('SH approximation of order 0')
sh0 = np.zeros(45)
sh0[0] = sh_signal[0]
signal_sphere = np.dot(sh0, B_highres)
print(sh0)
if interactive :
    my_visu(signal_sphere, sphere, False,
            title="Signal with SH order 0")

sh2 = np.zeros(45)
sh2[0:6] = sh_signal[0:6]
print(sh2)
print('SH approximation of order 2')
signal_sphere = np.dot(sh2, B_highres)
if interactive :
    my_visu(signal_sphere, sphere, False,
            title="Signal with SH order 2")


sh4 = np.zeros(45)
sh4[0:15] = sh_signal[0:15]
print(sh4)
print('SH approximation of order 4')
signal_sphere = np.dot(sh4, B_highres)
if interactive :
    my_visu(signal_sphere, sphere, False,
            title="Signal with SH order 4")

# alternatively, sh_to_sf and sf_to_sh functions are cool!
from dipy.reconst.shm import sh_to_sf, sf_to_sh
print(sh_signal)
sf_signal = sh_to_sf(sh_signal, sph_gtab, 8)
print(sf_signal)
print('SF reconstruction error:', np.mean(np.abs(signal_sph - sf_signal)))

sh_signal2 = sf_to_sh(sf_signal, sph_gtab, sh_order=8)
print(sh_signal2)
print('SH reconstruction error:', np.mean(np.abs(sh_signal - sh_signal2)))


#
# Part 2 of Demo - Playing with the diffusion ODF
#

# First, the single and multi-tensor models have an analytical diffusion ODF solution
print('Diffusion ODF')
odf = multi_tensor_odf(sphere.vertices, mevals, angles, fractions)

if interactive :
    my_visu(odf, sphere,
            title="Ground-truth population")

from dipy.reconst.shm import CsaOdfModel, QballModel
# qball model
qball_model = QballModel(gtab, 6)
fit = qball_model.fit(signal)
odf = fit.odf(sphere)
#if interactive :
#    my_visu(odf, sphere,
#            title="Q-Ball model Orientation Probabilities")

print('ODF min, max', odf.min(), odf.max())
odf = (odf - odf.min()) / (odf.max() - odf.min())
print('ODF_MinMaxNormalized min, max', odf.min(), odf.max())
if interactive :
    my_visu(odf, sphere,
            title="Q-Ball model Orientation Probabilities, normalized")

print('Qball GFA', fit.gfa)

# The q-ball diffusion ODF is not normalized and a smoothed version of the 'real' diffusion ODF

# CSA always best with sh_order 4
# otherwise, too noisy in the high frequencies, regularization would be needed
csa_model = CsaOdfModel(gtab, 4)
fit = csa_model.fit(signal)
odf = fit.odf(sphere)
if interactive :
    my_visu(odf, sphere,
            title="CSA model Orientation Probabilities")

odf = (odf - odf.min()) / (odf.max() - odf.min())
print('ODF_MinMaxNormalized min, max', odf.min(), odf.max())
if interactive :
    my_visu(odf, sphere,
            title="CSA model Orientation Probabilities, normalized")

print('CSA-Qball GFA', fit.gfa)
# The csa q-ball diffusion ODF is normalized and a closer approximation to the 'real' dODF
# However, vulnerable to noise and unstable for high SH order

# For more, see Dipy gallery
# http://nipy.org/dipy/examples_built/reconst_csa.html#example-reconst-csa

#
# Part 3 - Playing with Spherical Deconvolution 
#
print('Fiber ODFs with constrained spherical deconvolution (CSD)')
# First, we need a response function to deconvolve
# In real life, we estimate it from the data or we can fix it
# See http://nipy.org/dipy/examples_built/reconst_csd.html#example-reconst-csd

from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel, recursive_response

# Here, lets fix it and play with it.
response = (np.array([ 0.0015, 0.0004, 0.0004]), 1)
print(signal)
print(response)
# This is exactly the response that generated the synthetic signal. This is of course "cheating" a bit
# That is, a diffusion tensor [15,4,4] and mean normalized b0 signal 


csd_model = ConstrainedSphericalDeconvModel(gtab, response)
fit = csd_model.fit(signal)
fodf = fit.odf(sphere)

if interactive :
    my_visu(fodf, sphere,
            title="CSD Model fiber Orientation Probabilities")

from dipy.direction.peaks import peak_directions
directions, _, _ = peak_directions(fodf, sphere)
nufo = directions.shape[0]
print('Number of Fiber Orientations (NuFO):', nufo)

afd_total = fit.shm_coeff[0]
print('Total AFD:', afd_total)


# For more, see Dipy galleries
# http://nipy.org/dipy/examples_built/reconst_csd.html#example-reconst-csd
# http://nipy.org/dipy/examples_built/reconst_forecast.html#example-reconst-forecast

#
# Part 4 - Playing with multi-shell data
#
print('Multi-shell reconstruction')

grad = np.loadtxt('grad-multi-shell-100dirs.txt')
# Sorting gradients by shell, to ensure correct color assignment.
grad = grad[grad[:, 3].argsort()]
bvals = grad[:,3]
bvecs = grad[:,0:3]
print(bvals)

gtab = gradient_table(bvals, bvecs)

n_b0 = np.sum(bvals==0)
n_b300 = np.sum(bvals==300)
n_b1000 = np.sum(bvals==1000)
n_b2000 = np.sum(bvals==2000)
print('Number of b=0 images', n_b0)
print('Number of b=300 images', n_b300)
print('Number of b=1000 images', n_b1000)
print('Number of b=2000 images', n_b2000)

colors_b300 = window.colors.red * np.ones((n_b300,3))
colors_b1000 = window.colors.blue * np.ones((n_b1000,3))
colors_b2000 = window.colors.cyan * np.ones((n_b2000,3))
colors = np.vstack((colors_b300, colors_b1000, colors_b2000))
colors = np.ascontiguousarray(colors)

pts_actor = actor.point(gtab.gradients[~gtab.b0s_mask], colors, point_radius=100)
ren.add(pts_actor)
if interactive:
    window.show(ren, title="Multi-shell gradients distribution per shell",
                size=WINDOW_SIZE)
    ren.rm(pts_actor)

    ren.add(actor.point(gtab.bvecs[~gtab.b0s_mask], colors, point_radius=0.05))
    window.show(ren, title="Multi-shell gradients distribution, overall.",
                size=WINDOW_SIZE)


signal, sticks = multi_tensor(gtab, mevals, S0=S0, angles=angles,
                              fractions=fractions, snr=SNR)

from dipy.reconst.shore import ShoreModel

radial_order = 6
zeta = 700
lambdaN = 1e-8
lambdaL = 1e-8
asm = ShoreModel(gtab, radial_order=radial_order,
                 zeta=zeta, lambdaN=lambdaN, lambdaL=lambdaL)

asmfit = asm.fit(signal)

coeff = asmfit.shore_coeff
print(coeff.shape, coeff)

dodf = asmfit.odf(sphere)
if interactive :
    my_visu(dodf, sphere,
            title="SHORE model diffusion ODF")

print('Return to Origin Probability (RTOP):', asmfit.rtop_signal())
print('Mean Squared Displacement (MSD):', asmfit.msd())

# For more, see the Dipy galleries
# http://nipy.org/dipy/examples_built/reconst_csd.html#example-reconst-shore
# http://nipy.org/dipy/examples_built/reconst_forecast.html#example-reconst-shore-metrics
# http://nipy.org/dipy/examples_built/reconst_forecast.html#example-reconst-mapmri
 


# see multi-shell example in the mc-csd folder below
