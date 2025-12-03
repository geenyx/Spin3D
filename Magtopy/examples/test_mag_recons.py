#!/usr/bin/env python3
# Author: Marisel Di Pietro Martinez
# Any work making use of these codes should refer to, and cite, the
# following paper:
# Di Pietro Martínez et al., Phys. Rev. B 107, 094425 (2023)
# [https://doi.org/10.1103/PhysRevB.107.094425]

import numpy as np
import os
import matplotlib.pyplot as plt
from magtopy.magrecons import mag_recons_functions as mr

##############################################
############ PARAMETERS TO CHANGE ############
# Import data
FILE_IN = 'B1_projections.npz'
phis = np.load(FILE_IN)['phis']
thetas = np.load(FILE_IN, allow_pickle = True)['thetas']
proj = np.load(FILE_IN)['proj']
proj_mask = np.load(FILE_IN)['proj_mask'] #if there if one (optional)
L = proj.shape[0]

# Set output name
DIR_NAME = 'B1_recons'

# Import the original to compare (optional)
FILE_IN = 'B1.npz'
MX = np.load(FILE_IN)['MX']
MY = np.load(FILE_IN)['MY']
MZ = np.load(FILE_IN)['MZ']
ABS = np.load(FILE_IN)['ABS']
#MX = np.zeros((L,L,L))
#MY = np.zeros((L,L,L))
#MZ = np.zeros((L,L,L))
#ABS = np.zeros((L,L,L))

# Parameters for the reconstruction
total_ite = 20
n_ite = 2 #Print intermediate results each n_ite
step_x = 0.001
step_y = 0.001
step_z = 0.001
extras = {
    #"mask" : proj_mask # (optional)
    }

# Propose initial guess
# zeros
MX_guess = np.zeros((L,L,L), dtype=np.float32)
MY_guess = np.zeros((L,L,L), dtype=np.float32)
MZ_guess = np.zeros((L,L,L), dtype=np.float32)
##############################################
##############################################

# Convert angles to radians and define the rotation matrix elements
phis_r = np.float32(np.radians(phis))
thetas_r1 = np.float32(np.radians(thetas[0]))
thetas_r2 = np.float32(np.radians(thetas[1]))
thetas_r = np.array([thetas_r1, thetas_r2], dtype = object)
R31, R32, R33 = mr.rotations(phis_r, thetas_r)

# Create and change to new directory
DIR_NAME = mr.make_dir(DIR_NAME)
os.chdir(DIR_NAME)

##### Perform magnetic reconstruction in GPU #####
steps = []
error = []
for ite in range(total_ite//n_ite):
    MX_recons, MY_recons, MZ_recons, error_rel, last_steps, grad_x, grad_y, grad_z, proj_recons, step_x, step_y, step_z = mr.mag_recons(MX_guess, MY_guess, MZ_guess, phis_r, thetas_r, R31, R32, R33, ABS, proj, n_ite, step_x = step_x, step_y = step_y, step_z = step_z, **extras)
    
    MX_guess = MX_recons.copy()
    MY_guess = MY_recons.copy()
    MZ_guess = MZ_recons.copy()
    
    extras = {
        #"mask" : proj_mask, # (optional)
        "grad_x" : grad_x,
        "grad_y" : grad_y,
        "grad_z" : grad_z,
        "proj_guess" : proj_recons,
        "error_guess" : error_rel[-1]
    }
    
    error = np.concatenate((error, error_rel), axis=None)
    steps = np.concatenate((steps, last_steps), axis=None)

    ##############################################
    # Check some slices from reconstruction and compare to original
    slices = int(L*0.5)  
    plt.figure(figsize=(8, 6), dpi=100)
    plt.subplot(3,2,1)
    plt.imshow(MX[:,:,slices])
    plt.title('original')
    plt.ylabel('MX')
    plt.colorbar()
    plt.subplot(3,2,2)
    plt.imshow(MX_recons[:,:,slices])
    plt.title('guess')
    plt.colorbar()
    plt.subplot(3,2,3)
    plt.imshow(MY[:,:,slices])
    plt.ylabel('MY')
    plt.colorbar()
    plt.subplot(3,2,4)
    plt.imshow(MY_recons[:,:,slices])
    plt.colorbar()
    plt.subplot(3,2,5)
    plt.imshow(MZ[:,:,slices])
    plt.ylabel('MZ')
    plt.colorbar()
    plt.subplot(3,2,6)
    plt.imshow(MZ_recons[:,:,slices])
    plt.colorbar()
    plt.savefig("%s_magrecons_ite%d_xy.png" % (DIR_NAME, ite*n_ite), dpi=100, bbox_inches='tight')
    
    plt.figure(figsize=(8, 6), dpi=100)
    plt.subplot(3,2,1)
    plt.imshow(MX[:,slices,:])
    plt.title('original')
    plt.ylabel('MX')
    plt.colorbar()
    plt.subplot(3,2,2)
    plt.imshow(MX_recons[:,slices,:])
    plt.title('guess')
    plt.colorbar()
    plt.subplot(3,2,3)
    plt.imshow(MY[:,slices,:])
    plt.ylabel('MY')
    plt.colorbar()
    plt.subplot(3,2,4)
    plt.imshow(MY_recons[:,slices,:])
    plt.colorbar()
    plt.subplot(3,2,5)
    plt.imshow(MZ[:,slices,:])
    plt.ylabel('MZ')
    plt.colorbar()
    plt.subplot(3,2,6)
    plt.imshow(MZ_recons[:,slices,:])
    plt.colorbar()
    plt.savefig("%s_magrecons_ite%d_xz.png" % (DIR_NAME, ite*n_ite), dpi=100, bbox_inches='tight')
    
    plt.figure(figsize=(8, 6), dpi=100)
    plt.subplot(3,2,1)
    plt.imshow(MX[slices,:,:])
    plt.title('original')
    plt.ylabel('MX')
    plt.colorbar()
    plt.subplot(3,2,2)
    plt.imshow(MX_recons[slices,:,:])
    plt.title('guess')
    plt.colorbar()
    plt.subplot(3,2,3)
    plt.imshow(MY[slices,:,:])
    plt.ylabel('MY')
    plt.colorbar()
    plt.subplot(3,2,4)
    plt.imshow(MY_recons[slices,:,:])
    plt.colorbar()
    plt.subplot(3,2,5)
    plt.imshow(MZ[slices,:,:])
    plt.ylabel('MZ')
    plt.colorbar()
    plt.subplot(3,2,6)
    plt.imshow(MZ_recons[slices,:,:])
    plt.colorbar()
    plt.savefig("%s_magrecons_ite%d_yz.png" % (DIR_NAME, ite*n_ite), dpi=100, bbox_inches='tight')

    # Check the gradient
    plt.figure(figsize=(8, 6), dpi=100)
    plt.subplot(3,3,1)
    plt.imshow(grad_x[:,:,slices])
    plt.colorbar()
    plt.subplot(3,3,2)
    plt.imshow(grad_x[:,slices,:])
    plt.colorbar()
    plt.subplot(3,3,3)
    plt.imshow(grad_x[slices,:,:])
    plt.colorbar()
    plt.subplot(3,3,4)
    plt.imshow(grad_y[:,:,slices])
    plt.colorbar()
    plt.subplot(3,3,5)
    plt.imshow(grad_y[:,slices,:])
    plt.colorbar()
    plt.subplot(3,3,6)
    plt.imshow(grad_y[slices,:,:])
    plt.colorbar()
    plt.subplot(3,3,7)
    plt.imshow(grad_z[:,:,slices])
    plt.colorbar()
    plt.subplot(3,3,8)
    plt.imshow(grad_z[:,slices,:])
    plt.colorbar()
    plt.subplot(3,3,9)
    plt.imshow(grad_z[slices,:,:])
    plt.colorbar()

    plt.savefig("%s_magrecons_ite%d_grad.png" % (DIR_NAME, ite*n_ite), dpi=100, bbox_inches='tight')
    ##############################################
    ##############################################

# Save data into file
np.savez_compressed('%s_recons' % DIR_NAME, MX=MX_recons, MY=MY_recons, MZ=MZ_recons, error_rel=error, steps=steps, grad_x=grad_x, grad_y=grad_y, grad_z=grad_z, proj_recons = proj_recons)

# Check error reduction
plt.figure(figsize=(8, 6), dpi=100)
plt.plot(error, 'ro-', linewidth=1.5)
plt.ylabel('error [%]')
plt.xlabel('iterations')
plt.ylim([0,100])
plt.savefig("%s_magrecons_error.png" % DIR_NAME, dpi=100, bbox_inches='tight')

# See step variation
plt.figure(figsize=(8, 6), dpi=100)
plt.plot(steps, 'ro-', linewidth=1.5)
plt.ylabel('step_x')
plt.xlabel('iterations')
plt.savefig("%s_magrecons_steps.png" % DIR_NAME, dpi=100, bbox_inches='tight')
