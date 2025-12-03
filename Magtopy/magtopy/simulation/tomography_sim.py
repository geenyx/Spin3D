#!/usr/bin/env python3
import numpy as np
from scipy.fft import fft2
from numpy.fft import fftshift

def cdi_intensity(proj, DynRange = 0, Poisson = False):
    """
    Calculates the intensities that would be measured in a CDI experiment.
    
    Parameters
    ----------
    proj : ndarray
        Array of projections (L,L,n_angles)
    DynRange : float32
        Dynamic range. Default is 0.
    Poisson : bool
        Whether the measurement has (Poisson) noise or not. Default is False.
    
    Returns
    -------
    I : ndarray (int)
        Array with the intensities (counts) measured with CDI (L,L,n_angles)
    """

    FT = fftshift(fft2(proj, axes=(0,1)), axes=(0,1))
    I = (np.real(FT))**2+(np.imag(FT))**2
    
    if DynRange > 0:
        I_measured = (DynRange * (I/np.max(I))).astype(int)
    else:
        I_measured = I.astype(int)
    if Poisson:
        I_measured = np.random.poisson(I_measured)
    
    return I_measured
