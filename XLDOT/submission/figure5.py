import magpack
from ThreeDViewer import image
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import functions
import numpy as np


def plot_histograms():
    for regime in ['smoothly_varying']:
        vf = functions.load_reference(regime)
        mask = magpack.vectorop.magnitude(vf) > 0
        div = np.abs(magpack.vectorop.divergence(vf))[mask]
        curl = magpack.vectorop.magnitude(magpack.vectorop.curl(vf))[mask]
        for method in ['dual_axis', 'single_axis']:
            angular_error = functions.load_errors(regime, method)[0][mask]
            n_bins = 9
            stop = 0.8
            data = np.zeros((n_bins, n_bins))
            equiv_hist = np.zeros((n_bins, n_bins))
            intervals = np.linspace(0, stop, n_bins)
            for ii, (x_edge_low, x_edge_high) in enumerate(zip(intervals, intervals + 0.1)):
                div_mask = np.logical_and(div > x_edge_low, div <= x_edge_high)
                for jj, (y_edge_low, y_edge_high) in enumerate(zip(intervals, intervals + 0.1)):
                    curl_mask = np.logical_and(curl > y_edge_low, curl <= y_edge_high)
                    intersection = np.logical_and(div_mask, curl_mask)
                    equiv_hist[ii, jj] = np.count_nonzero(intersection)

                    data[ii, jj] = np.nanmedian(angular_error[intersection])
            plt.figure()
            im = plt.imshow(equiv_hist.T, interpolation='nearest', origin='lower',
                            extent=(0., stop + 0.1, 0., stop + 0.1), norm=LogNorm())
            plt.ylabel('curl')
            plt.xlabel('divergence')
            plt.title("Bin count with specified curl and divergence")
            plt.colorbar(im)
            plt.show()

            plt.figure()
            plt.imshow(data.T, interpolation='nearest', origin='lower', extent=(0., stop + 0.1, 0., stop + 0.1), vmin=0,
                       vmax=45)
            plt.ylabel('curl')
            plt.xlabel('divergence')
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Average angular error / °')
            plt.show()


def plot_examples():
    methods = ['single_axis', 'triple_axis']
    regime = 'smoothly_varying'
    roi = (..., slice(95, 115), slice(94, 114), 113)

    # the feature is at an angle to the conventional slicing axes so needs to be rotated
    ref_vf = functions.load_reference(regime)
    ref_vf_rot = magpack.rotations.rotate_vector_field(ref_vf, magpack.rotations.roty(-30), order=1)
    curl_rot = magpack.vectorop.magnitude(magpack.vectorop.curl(ref_vf_rot))
    div_rot = np.abs(magpack.vectorop.divergence(ref_vf_rot))

    plt.imshow(div_rot[roi].T, origin='lower')
    plt.title("Divergence")
    plt.colorbar()
    plt.show()

    plt.imshow(curl_rot[roi].T, origin='lower')
    plt.title("Curl")
    plt.colorbar()
    plt.show()

    for method in methods:
        vf = functions.load_recons(regime, method)
        vf_rot = magpack.rotations.rotate_vector_field(vf, magpack.rotations.roty(-30), order=1)

        image.plot_quiver(vf_rot[roi], axial=True, skip=1, mode=1, title=f'Reconstruction ({method})')
        image.plot_quiver(ref_vf_rot[roi], axial=True, skip=1, mode=1, title=f'Reference')


if __name__ == '__main__':
    plot_histograms()
    plot_examples()
