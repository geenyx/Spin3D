from itertools import product

import magpack.vectorop
import numpy as np
from matplotlib import pyplot as plt
import functions


def averaging_trend():
    regime = 'smoothly_varying'
    vf_ref = functions.load_reference(regime)
    mask = magpack.vectorop.magnitude(vf_ref) > 0

    x_axis = np.array([5, 10, 15, 20])
    data = np.empty((x_axis.shape[0], 5))
    for ii, avg_number in enumerate([5, 10, 15, 20]):
        reconstruction = magpack.io.load_mat(f'output/{regime}/triple_axis/average{avg_number:02d}.mat')
        uncertainty = reconstruction['uncertainty']
        data[ii, :] = np.percentile(uncertainty[mask], [5, 25, 50, 75, 95])

    fig, ax = plt.subplots()
    ax.fill_between(x_axis, data[:, 0], data[:, 4], color=(0.75, 0.75, 0.75, 1))
    ax.fill_between(x_axis, data[:, 1], data[:, 3], color=(0.5, 0.5, 0.5, 1))

    ax.plot(x_axis, data[:, 2], color='k', label="median", marker='o')
    ax.margins(x=0)
    plt.show()


def uncertainty_plots():
    roi_smooth = (..., slice(50, 80), slice(50, 90), 74)
    roi_granular = (..., slice(65, 95), slice(55, 95), 54)
    regimes = ['smoothly_varying', 'granular']
    methods = ['single_axis', 'dual_axis', 'triple_axis']

    for (regime, roi), method in product(zip(regimes, [roi_smooth, roi_granular]), methods):
        angular_error, uncertainty = functions.load_errors(regime, method)
        plt.imshow(angular_error[25:-25, 25:-25, roi[-1]].T, vmin=0, vmax=90, cmap='turbo', origin='lower')
        plt.title(f"Angular error {regime} ({method})")
        plt.colorbar()
        plt.show()

        plt.imshow(uncertainty[25:-25, 25:-25, roi[-1]].T, vmin=0, vmax=60, cmap='afmhot', origin='lower')
        plt.title(f"Uncertainty {regime} ({method})")
        plt.colorbar()
        plt.show()


if __name__ == '__main__':
    averaging_trend()
    uncertainty_plots()