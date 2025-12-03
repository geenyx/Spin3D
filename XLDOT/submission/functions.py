import magpack.io as io
import numpy as np


def _load_reconstruction_data(regime, method):
    if regime not in ['smoothly_varying', 'granular']:
        raise ValueError("Regime must be either 'smoothly_varying' or 'granular'.")
    if method not in ['single_axis', 'dual_axis', 'triple_axis']:
        raise ValueError("Method must be either 'single_axis', 'dual_axis' or 'triple_axis'.")
    return io.load_mat(f'output/{regime}/{method}/average20.mat')


def load_recons(regime, method):
    data = _load_reconstruction_data(regime, method)
    return np.stack([data['final_x'], data['final_y'], data['final_z']], axis=0)


def load_errors(regime, method):
    data = _load_reconstruction_data(regime, method)
    return data['theta_error'], data['uncertainty']


def load_reference(regime):
    if regime not in ['smoothly_varying', 'granular']:
        raise ValueError("Regime must be either 'smoothly_varying' or 'granular'.")
    data = io.load_mat(f'reference_structures/{regime}.mat')
    return np.stack([data['mx'], data['my'], data['mz']], axis=0)

