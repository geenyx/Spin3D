#!/usr/bin/env python3
# Creates a 3D magnetic sample
import numpy as np
from math import cos, sin

def make_sample(L, shape='cylinder'):
    """
    Creates a LxLxL support with a sample.
    
    Parameters
    ----------
    L : int
        Size of the support LxLxL
    shape : 'cylinder' (default), 'cube'
        The shape of the sample
    
    Return
    ------
    ABS : ndarray
        Mask with 0 where there's no sample and 1 where there is.
    """
    ABS = np.zeros((L,L,L), dtype=np.float32)

    center = int(L/2)
    r_cut = int(L/4)
    r_cut2 = r_cut*r_cut
    height = int(L/4)
    
    if shape=='cylinder':
        for x in range(L):
            xc = x - center
            for y in range(L):
                yc = y - center
                for z in range(L):
                    zc = z - center
                    rr = int(xc*xc+yc*yc)
                    if (rr < r_cut2 and abs(zc) < height):
                        ABS[x,y,z] = 1
    elif shape=='cube':
        ABS[r_cut:-r_cut,r_cut:-r_cut,r_cut:-r_cut] = 1
    
    return ABS

def add_charge(ABS, F_CHARGE=1, method='constant'):
    """
    Adds electronic charge.
    
    Parameters
    ----------
    ABS : ndarray
        It is equal to 1 where there's sample and zero where there's not
    
    Return
    ------
    RHO : ndarray
        Electronic density
    """
    
    RHO = np.zeros_like(ABS)
    if method == 'constant':
        RHO = F_CHARGE*ABS
    
    return RHO

def add_mag(ABS, F_MAG=1, method='vortex', lamb=5, phi=0, theta=0):
    """
    Adds 3D magnetization to the sample
    
    Parameters
    ----------
    ABS : ndarray
        It is equal to 1 where there's sample and zero where there's not
    lamb : float32 (optional, for method=stripes)
        Wave length
    phi, theta : float32 (optional, for method=stripes)
        Angles in degrees that indicate the direction of the spins
    
    Return
    ------
    MX, MY, MZ: ndarrays
        The three components of the magnetization
    """
    
    MX = np.zeros_like(ABS)
    MY = np.zeros_like(ABS)
    MZ = np.zeros_like(ABS)
    
    L = ABS.shape[0]

    center = int(L/2)
    r_cut2 = int(L*L/16)
    height = int(L/4)
    
    if method=='vortex':
        for x in range(L):
            xc = x - center
            for y in range(L):
                yc = y - center
                for z in range(L):
                    zc = z - center
                    MX[x,y,z] = -np.sign(yc)
                    MY[x,y,z] = np.sign(xc)
    elif method=='stripes':
        phi_r = np.radians(phi)
        theta_r = np.radians(theta)
        RX = -sin(phi_r)*cos(theta_r)
        RY = -sin(phi_r)*sin(theta_r)
        RZ = cos(phi_r)

        for x in range(L):
            xc = x - center
            stripe = 2 * (int(x/lamb)%2) - 1
            stripe_x = stripe*RX
            stripe_y = stripe*RY
            stripe_z = stripe*RZ
            for y in range(L):
                yc = y - center
                for z in range(L):
                    zc = z - center
                    r = int(xc*xc+yc*yc)
                    if (r < r_cut2 and abs(zc) < height):
                        ABS[x,y,z] = 1
                        MX[x,y,z] = stripe_x
                        MY[x,y,z] = stripe_y
                        MZ[x,y,z] = stripe_z
    #elif method=='random':
        #TODO
        # = np.random.randint(2, size=10)
    
    MX = F_MAG*MX*ABS
    MY = F_MAG*MY*ABS
    MZ = F_MAG*MZ*ABS
    
    return MX, MY, MZ
