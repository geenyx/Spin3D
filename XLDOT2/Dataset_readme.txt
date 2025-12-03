2024 Dec 02
X-ray linear dichroic orientation tomography data and reconstruction codes

On-resonance datasets:
data_lh_shift: projections obtained using linear horizontal (LH) polarisation, shifted for 3D alignment
data_lv_shift: projections obtained using linear vertical (LV) polarisation, shifted for 3D alignment
theta_lh: rotation angle at which the LH projection was measured
theta_lv: rotation angle at which the LV projection was measured
R_lh: rotation matrix (accounting for sample tilting) describing the sample orientation for LH projection 
R_lv: rotation matrix (accounting for sample tilting) describing the sample orientation for LV projection 

aligned_proj_tilt1.mat (data from reference tilt, 0 degree tilt)
-> aligned_r: [structure]
	-> data_lh_shift: [760×650×280 single]
        -> data_lv_shift: [760×650×280 single]
        -> R_lh: [3×3×280 double]
        -> R_lv: [3×3×280 double]
        -> theta_lh: [1×280 single]
        -> theta_lv: [1×280 single]

aligned_proj_tilt2.mat (data from first tilt, 30 degree tilt along z)
-> aligned_a: [structure]
	-> data_lh_shift: [760×650×281 single]
        -> data_lv_shift: [760×650×281 single]
        -> R_lh: [3×3×281 double]
        -> R_lv: [3×3×281 double]
        -> theta_lh: [1×281 single]
        -> theta_lv: [1×281 single]

aligned_proj_tilt3.mat (data from second tilt, 30 degree tilt along x)
-> aligned_c: [structure]
	-> data_lh_shift: [760×650×281 single]
        -> data_lv_shift: [760×650×281 single]
        -> R_lh: [3×3×281 double]
        -> R_lv: [3×3×281 double]
        -> theta_lh: [1×281 single]
        -> theta_lv: [1×281 single]

aligned_proj_tilt4.mat (data from third tilt, -30 degree tilt along x)
-> aligned_c: [structure]
	-> data_lh_shift: [760×650×281 single]
        -> data_lv_shift: [760×650×281 single]
        -> R_lh: [3×3×281 double]
        -> R_lv: [3×3×281 double]
        -> theta_lh: [1×281 single]
        -> theta_lv: [1×281 single]


Off-resonance dataset:
data_ch_shift: tomographic projections obtained off-resonance, shifted for 3D alignment
theta_ch: rotation angle at which the projection was measured
R_ch: rotatin matrix (accounting for sample tilting) describing the sample orientation
factor_edensity: conversion factor (from phase to electron density)

off_resonance_projections.mat
-> aligned_q: [structure]
	-> data_ch_shift: [760×650×350 double]
        -> R_ch: [3×3×350 double]
        -> theta_ch: [1×350 single]
-> factor_edensity: [double]


Reconstructions:
electron_density_tomogram.mat
electron_density: the electron density reconstruction in units of number of electrons per cubic Angstrom
-> electron_density: [650x650x760 single]

xldot-reconstruction.mat
-> mean_x: [330x330x540 single] (x component of the orientation)
-> mean_y: [330x330x540 single] (y component of the orientation)
-> mean_z: [330x330x540 single] (z component of the orientation)

reconstruction.vtk
-> scalar fields: 
	-> eden (electron density)
	-> segmentation (segmentation index)
-> vector fields:
	-> mag (orientation of the c-axis)


Statistics:
segmentation_stats.csv (contains grain segmentation statistics such as volume, surface area, sphericity, etc. First row reserved for headers describing each column)


Scripts:
reconstruction.m: 
(Section 1.0) Loads projections
(Section 2.0) Performs scalar reconstruction to obtain electron density tomogram and defines a sample mask
(Section 3.1) Computes the difference of projections (LH - LV) and combines them to a single stack
(Section 3.2) Perform reconstructions
(Section 3.3) Average reconstructions


