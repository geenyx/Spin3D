#!/usr/bin/env python3
# Author: Marisel Di Pietro Martinez
# Any work making use of these codes should refer to, and cite, the
# following paper:
# Di Pietro Martínez et al., Phys. Rev. B 107, 094425 (2023)
# [https://doi.org/10.1103/PhysRevB.107.094425]

# third party packages
import numpy as np
import math
import os

# third party packages CUDA
import pycuda.autoinit
import pycuda.gpuarray as gpuarray

# local packages
from magtopy.magrecons import mag_recons_pycuda_lib as mrc

def root_mean_square(XGR, X):
    # Normalized root mean square error: NRMSE
    L = X.shape[0]
    Linv = 1.0/L
    NRMSE = np.zeros(L)
    for z in range(L):
        xmax = np.max(X[:,:,z])
        xmin = np.min(X[:,:,z])
        xdif = 1.0/(xmax - xmin)
        NRMSE[z] = xdif * np.sqrt(Linv * Linv * np.sum(np.square(XGR[:,:,z] - X[:,:,z])))
    return NRMSE

def root_mean_square_ite(XGR, X):
    # Normalized root mean square error: NRMSE
    L = X.shape[0]
    Linv = 1.0/L
    xmax = np.max(X)
    xmin = np.min(X)
    xdif = 1.0/(xmax - xmin)
    NRMSE = xdif * np.sqrt(Linv * Linv * Linv * np.sum(np.square(XGR - X)))
    return NRMSE

def make_dir(FILE_NAME):
    """
    Handle making a new directory for the data.
    
    :FILE_NAME: string with the desired file name
    """
    fileok = False
    while not fileok:
        try:
            os.mkdir(FILE_NAME)
        except:
            x = input('File %s already exist. Rename? [(y)/n/c] (n: overwrite, c: cancel)' % FILE_NAME)
            if x=='c':
                print('Quiting.')
                quit()
            elif x=='n':
                print('Overwriting.')
                fileok = True
            else:
                FILE_NAME = input('File Name: ')
        else:
            fileok = True
    return FILE_NAME

def rotations(phis_rad, thetas_rad):
    """
    Parameters
    ----------
    phis_rad : ndarray
        Angles in radians to rotate clockwise around z axis
    thetas_rad : ndarray
        Angles in radians to rotate counter-clockwise around y axis
    
    Returns
    -------
    R31, R32, R33 : ndarrays
        The elements R31, R32, R33 of the rotation matrix R.
    """
    n_angles = 0
    for i in range(phis_rad.size):
        n_angles += thetas_rad[i].size
    R31 = np.zeros(n_angles, dtype=np.float32)
    R32 = np.zeros(n_angles, dtype=np.float32)
    R33 = np.zeros(n_angles, dtype=np.float32)
    angle = 0
    angleset = 0
    for phi in phis_rad:
        cosphi = round(math.cos(phi),5)
        sinphi = round(math.sin(phi),5)
        for theta in thetas_rad[angleset]:
            costheta = round(math.cos(theta),5)
            sintheta = round(math.sin(theta),5)
            R31[angle] = - cosphi*sintheta
            R32[angle] = - sinphi*sintheta
            R33[angle] = costheta
            angle += 1
        angleset += 1
    return R31, R32, R33

def project(MX_d, MY_d, MZ_d, RHO_d, phis_rad, thetas_rad, L, n_angles, R3_1, R3_2, R3_3, mask_d):
    """
    Compute magnetic projections along z-axis
    
    Parameters
    ----------
    MX_d, MY_d, MZ_d : gpuarrays
        Magnetization components
    RHO_d : gpuarray
        Electronic charge
    phis_rad : ndarray
        Angles in radians to rotate clockwise around z axis
    thetas_rad : ndarray
        Angles in radians to rotate counter-clockwise around y axis
    L : int
        Size of the cube LxLxL
    n_angles : int
        Number of angles in total (phis_rad.size*thetas_rad.size)
    R3_1, R3_2, R3_3 : ndarray
        The elements R31, R32, R33 of the rotation matrix R.
    mask_d : gpuarrays
        Mask for the projections
        
    Returns
    -------
    proj_d : gpuarray
        Magnetic projection along z for each rotation of the sample
    """
    proj_d = gpuarray.zeros((L,L,n_angles), dtype=np.float32)
    angle = 0
    angleset = 0
    for phi in phis_rad:
        for theta in thetas_rad[angleset]:
            R31 = R3_1.flat[angle]
            R32 = R3_2.flat[angle]
            R33 = R3_3.flat[angle]
            MZ_rot_d = mrc.mag_rotate(MX_d, MY_d, MZ_d, RHO_d, R31, R32, R33)
            MZ_rot_d = mrc.rot_z(MZ_rot_d, phi)
            MZ_rot_d = mrc.rot_y(MZ_rot_d, theta)
            one_proj_d = mrc.project(MZ_rot_d)
            proj_d[:,:,angle] = one_proj_d.get()
            angle += 1
        angleset += 1
    return proj_d*mask_d

def gradient(dif_proj_d, phis_rad, thetas_rad, L, R3_1, R3_2, R3_3, ABS_d):
    """
    Compute gradient to perform gradient descent algorithm
    
    Parameters
    ----------
    dif_proj_d : gpuarray
        Difference between the guess and the measured projections
    phis_rad : ndarray
        Angles in radians to rotate clockwise around z axis
    thetas_rad : ndarray
        Angles in radians to rotate counter-clockwise around y axis
    L : int
        Size of the cube LxLxL
    R3_1, R3_2, R3_3 : ndarrays
        The elements R31, R32, R33 of the rotation matrix R
    ABS_d : gpuarray 
        Sample support
    
    Returns
    -------
    grad_x_d, grad_y_d, grad_z_d: gpuarrays
        Gradient to direct the update of the magnetization guess
    """
    grad_x_d = gpuarray.zeros((L,L,L), dtype=np.float32)
    grad_y_d = gpuarray.zeros((L,L,L), dtype=np.float32)
    grad_z_d = gpuarray.zeros((L,L,L), dtype=np.float32)
    angle = 0
    angleset = 0
    #dif_proj_h = dif_proj_d.get()
    for phi in phis_rad:
        for theta in thetas_rad[angleset]:
            # P^-1
            back_proj_d = mrc.back_project((dif_proj_d[:,:,angle]).copy())
            
            # R^-1
            back_proj_d = mrc.rot_y(back_proj_d, -theta)
            back_proj_d = mrc.rot_z(back_proj_d, -phi)
            
            R31 = R3_1.flat[angle]
            R32 = R3_2.flat[angle]
            R33 = R3_3.flat[angle]
            grad_x_d = grad_x_d.mul_add(1.0, back_proj_d, R31)
            grad_y_d = grad_y_d.mul_add(1.0, back_proj_d, R32)
            grad_z_d = grad_z_d.mul_add(1.0, back_proj_d, R33)
            
            angle += 1
        angleset += 1
            
    grad_x_d = 2*grad_x_d*ABS_d
    grad_y_d = 2*grad_y_d*ABS_d
    grad_z_d = 2*grad_z_d*ABS_d
    
    return grad_x_d, grad_y_d, grad_z_d

def mag_recons(MX_guess, MY_guess, MZ_guess, phis_rad, thetas_rad, R3_1, R3_2, R3_3, ABS, proj, n_ite = 10, step_x = 0.0001, step_y = 0.0001, step_z = 0.0001, max_ite = 10, val = 10, **kwargs):
    """
    Reconstruct the three components of the magnetization in a sample using the tomographic projections.
    
    Parameters
    ----------
    MX_guess, MY_guess, MZ_guess : ndarrays
        Initial guess for the magnetization components
    phis_rad : ndarray
        Angles in radians to rotate clockwise around z axis
    thetas_rad : ndarray
        Angles in radians to rotate counter-clockwise around y axis
    R3_1, R3_2, R3_3 : ndarrays
        The elements R31, R32, R33 of the rotation matrix R
    ABS : ndarray 
        Sample support
    proj : ndarray
        Measured tomographic projections
    n_ite : int
        Number of iterations to perform
    step : float32
        Step for the gradient descent.

    Returns
    -------
    MX_recons, MY_recons, MZ_recons : ndarrays
        The three components of the reconstructed magnetization
    error_rel : ndarray
        The percentual error vs iterations
    grad_x, grad_y, grad_z : ndarrays
        Last gradient calculated
    """

    # Copy to GPU
    ABS_d = gpuarray.to_gpu(np.float32(ABS))
    proj_d = gpuarray.to_gpu(np.float32(proj))
    MX_guess_d = gpuarray.to_gpu(MX_guess)
    MY_guess_d = gpuarray.to_gpu(MY_guess)
    MZ_guess_d = gpuarray.to_gpu(MZ_guess)
    #TODO agregar la reconstruction de rho
    RHO_d = gpuarray.zeros_like(MX_guess_d)

    sum_proj_d = gpuarray.sum(proj_d*proj_d, dtype=np.float32)
    sum_proj = sum_proj_d.get()
    
    L = proj.shape[0]
    n_angles = 0
    for i in range(phis_rad.size):
        n_angles += thetas_rad[i].size
    error_rel = np.zeros(n_ite)
    steps = np.zeros(n_ite)
    
    ##############################
    # Check extra arguments
    if ("mask" in kwargs):
        mask = kwargs["mask"]
        mask_d = gpuarray.to_gpu(np.float32(mask))
    else:
        mask = np.ones((L,L,n_angles))
        mask_d = gpuarray.to_gpu(np.float32(mask))
        
    if "proj_guess" in kwargs:
        proj_guess = kwargs["proj_guess"]
        proj_guess_d = gpuarray.to_gpu(np.float32(proj_guess))
    else:    
        proj_guess_d = project(MX_guess_d, MY_guess_d, MZ_guess_d, RHO_d, phis_rad, thetas_rad, L, n_angles, R3_1, R3_2, R3_3, mask_d)
        print('Calculate proj...')
    
    if "error_guess" in kwargs:
        error_guess = kwargs["error_guess"]
    else:
        dif_proj_d = proj_guess_d - proj_d
        dif_proj2_d = dif_proj_d*dif_proj_d
        error_guess_d = gpuarray.sum(dif_proj2_d, dtype=np.float32)
        error_guess = error_guess_d.get()*100/sum_proj
        print('Calculate error...')
    
    if ("grad_x" in kwargs) and ("grad_y" in kwargs) and ("grad_z" in kwargs):
        grad_x = kwargs["grad_x"]
        grad_x_d = gpuarray.to_gpu(np.float32(grad_x))
        grad_y = kwargs["grad_y"]
        grad_y_d = gpuarray.to_gpu(np.float32(grad_y))
        grad_z = kwargs["grad_z"]
        grad_z_d = gpuarray.to_gpu(np.float32(grad_z))
    else:
        grad_x_d, grad_y_d, grad_z_d = gradient(dif_proj_d, phis_rad, thetas_rad, L, R3_1, R3_2, R3_3, ABS_d)
        print('Calculate grad...')

    ##############################
    # Magnetic Reconstruction    
    for ite in range(n_ite):
        ite_step = 0
        label = 'check'
        while(ite_step < max_ite):
            # Gradient descent
            MX_new_d = MX_guess_d.mul_add(1.0, grad_x_d, np.float32(-step_x))
            MY_new_d = MY_guess_d.mul_add(1.0, grad_y_d, np.float32(-step_y))
            MZ_new_d = MZ_guess_d.mul_add(1.0, grad_z_d, np.float32(-step_z))

            ########## Calculate projections from the candidate###########
            proj_new_d = project(MX_new_d, MY_new_d, MZ_new_d, RHO_d, phis_rad, thetas_rad, L, n_angles, R3_1, R3_2, R3_3, mask_d)
            ###########################################################

            ########## Calculate error ################################
            dif_proj_d = proj_new_d - proj_d
            dif_proj2_d = dif_proj_d*dif_proj_d
            error_new_d = gpuarray.sum(dif_proj2_d, dtype=np.float32)
            error_new = error_new_d.get()*100/sum_proj
            ###########################################################

            ########## Test chosen step ###############################
            # Compare error
            #print(error_guess, error_new, sum_proj, error_new_d.get(), MZ_new_d.get()[57,57,57], proj_guess_d.get()[57,57,0], proj_new_d.get()[57,57,0], proj[57,57,0])
            derror = error_new - error_guess
            # kind of the Armijo rule...
            print('derror = %lf' % derror)
            if derror > 0:
                # too big!
                step_x = step_x/2
                step_y = step_y/2
                step_z = step_z/2
                ite_step += 1
                print('The step is too big!')
                label = 'less'
            elif abs(derror) != 0 and abs(derror) < val and ite_step < (max_ite - 1) and label != 'less':
                # it could be bigger...
                step_x = step_x*2
                step_y = step_y*2
                step_z = step_z*2
                ite_step += 1
                print('The step could be bigger...')
                label = 'more'
            else:
                if abs(derror) == 0:
                    print('derror is 0.')
                # The step is good!
                MX_guess_d = MX_new_d.copy()
                MY_guess_d = MY_new_d.copy()
                MZ_guess_d = MZ_new_d.copy()
                proj_guess_d = proj_new_d.copy()
                error_guess = error_new.copy()
                ite_step = max_ite
                print('The step is OK.')
                label = 'The step is OK.'
                #break

        ########## Calculate error gradient from the guess#########
        grad_x_d, grad_y_d, grad_z_d = gradient(dif_proj_d, phis_rad, thetas_rad, L, R3_1, R3_2, R3_3, ABS_d)
        ###########################################################
        
        error_rel[ite] = error_guess.copy()
        steps[ite] = step_x
        print('%d %lf %lf %d\n' % (ite, error_rel[ite], step_x, ite_step))

    # Copy to CPU
    MX_recons = MX_guess_d.get()
    MY_recons = MY_guess_d.get()
    MZ_recons = MZ_guess_d.get()
    grad_x = grad_x_d.get()
    grad_y = grad_y_d.get()
    grad_z = grad_z_d.get()
    proj_recons = proj_guess_d.get()
    
    return MX_recons, MY_recons, MZ_recons, error_rel, steps, grad_x, grad_y, grad_z, proj_recons, step_x, step_y, step_z
