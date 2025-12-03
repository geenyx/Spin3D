import functions
from ThreeDViewer import image
from itertools import product
from matplotlib import pyplot as plt

if __name__ == '__main__':
    roi_smooth = (..., slice(50, 80), slice(50, 90), 74)
    roi_granular = (..., slice(65, 95), slice(55, 95), 54)
    regimes = ['smoothly_varying', 'granular']
    methods = ['single_axis', 'dual_axis', 'triple_axis']

    for (regime, roi), method in product(zip(regimes, [roi_smooth, roi_granular]), methods):
        ref_vf = functions.load_reference(regime)
        vf = functions.load_recons(regime, method)
        angular_error = functions.load_errors(regime, method)[0]
        image.plot_quiver(ref_vf[..., roi[-1]], slice_axis=2, axial=True, mode=1,
                          title=f'Reference {regime}')
        image.plot_quiver(vf[..., roi[-1]], slice_axis=2, axial=True, mode=1,
                          title=f'Reconstruction {regime} ({method})')

        plt.imshow(angular_error[25:-25, 25:-25, roi[-1]].T, vmin=0, vmax=90, cmap='turbo', origin='lower')
        plt.title(f"Angular error {regime} ({method})")
        plt.colorbar()
        plt.show()

        image.plot_quiver(ref_vf[roi], slice_axis=2, axial=True, mode=1, skip=1, title=f'Reference {regime}')
        image.plot_quiver(vf[roi], slice_axis=2, axial=True, mode=1, skip=1,
                          title=f'Reconstruction {regime} ({method})')
