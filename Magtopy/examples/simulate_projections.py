#!/usr/bin/env python3
# Author: Marisel Di Pietro Martinez
# Any work making use of these codes should refer to, and cite, the
# following paper:
# Di Pietro Martínez et al., Phys. Rev. B 107, 094425 (2023)
# [https://doi.org/10.1103/PhysRevB.107.094425]

import pycuda.gpuarray as gpuarray
import numpy as np
from magtopy.magrecons import mag_recons_functions as mr

# Read sample data
FILE_NAME = 'B1.npz'
MX = np.load(FILE_NAME)['MX']
MY = np.load(FILE_NAME)['MY']
MZ = np.load(FILE_NAME)['MZ']
ABS = np.load(FILE_NAME)['ABS']
RHO = np.load(FILE_NAME)['RHO']
L = MX.shape[0]

# Define projection angles in degrees
phis = np.array([0, 90])
thetas1 = np.arange(-90,91,10)
thetas = np.array([thetas1, thetas1])
print(phis, thetas)

# Convert angles to radians and define the rotation matrix elements
phis_rad = np.float32(np.radians(phis))
thetas_rad = np.float32(np.radians(thetas))
R31, R32, R33 = mr.rotations(phis_rad, thetas_rad)

# Prepare data to pass to CUDA function
MX_d = gpuarray.to_gpu(np.float32(MX))
MY_d = gpuarray.to_gpu(np.float32(MY))
MZ_d = gpuarray.to_gpu(np.float32(MZ))
RHO_d = gpuarray.to_gpu(np.float32(RHO))
n_angles = thetas[0].size + thetas[1].size

# Load a non-trivial mask if required
mask = np.ones((L,L,n_angles))
mask_d = gpuarray.to_gpu(np.float32(mask))
print(L, mask.shape)

##### Calculate projections in GPU #####
proj_d = mr.project(MX_d, MY_d, MZ_d, RHO_d, phis_rad, thetas_rad, L, n_angles, R31, R32, R33, mask_d)

# Copy data to CPU
proj = proj_d.get()

np.savez_compressed('B1_projections', phis=phis, thetas=thetas, proj=proj, proj_mask = mask)
