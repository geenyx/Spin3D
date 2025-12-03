import functions
from ThreeDViewer import image

if __name__ == '__main__':
    roi_smooth = (slice(None, None), slice(50, 80), slice(50, 90), slice(74, 74 + 1))
    roi_granular = (slice(None, None), slice(65, 95), slice(55, 95), slice(54, 54 + 1))

    for regime, slice_idx, roi in zip(['smoothly_varying', 'granular'], [74, 54], [roi_smooth, roi_granular]):
        ref_vf = functions.load_reference(regime)
        vf = functions.load_recons(regime, 'triple_axis')
        image.plot_quiver(ref_vf[..., slice_idx], slice_axis=2, axial=True, mode=1, title='Reference')
        image.plot_quiver(vf[..., slice_idx], slice_axis=2, axial=True, mode=1, title='Reconstruction')

        image.plot_quiver(ref_vf[roi].squeeze(), slice_axis=2, skip=1, axial=True, mode=1, title='Reference')
        image.plot_quiver(vf[roi].squeeze(), slice_axis=2, skip=1, axial=True, mode=1, title='Reconstruction')
