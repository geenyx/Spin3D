An instruction manual for magnetic tomography reconstruction algorithms.

These codes can be used for the reconstruction of a magnetisation vector 
field from a set of single polarisation XMCD projections.
Each projection is defined by a single rotation matrix, which defines the 
orientation of the structure with respect to the X-ray beam. In this way, 
the code can deal with "arbitrary projections" and thus is not limited to a 
tomographic geometry.

Any work making use of these codes should cite the following paper:

Tomographic reconstruction of a three-dimensional magnetization vector 
field
Claire Donnelly, Sebastian Gliga, Valerio Scagnoli, Mirko Holler, J ̈org 
Raabe, Laura J. Heyderman and Manuel Guizar-Sicairos.
New Journal of Physics (2018)

The paper provides a description of the mathematical background of the 
codes, as well as a detailed analysis of the reconstruction of the provided 
micromagnetic simulation.


Quick overview of the codes provided:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
An example micromagnetic simulation of a pillar is provided: 
    cylinder_end.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XMCD projections (single circular left polarisation assumed) can be 
calculated from this simulation (or your own simulation) using the script:
    mag_tomo_calculate_projections.m
In this script, the projections are defined by a set of rotation matrices. 
An example calculation of rotation matrices is given for the dual-axis 
geometry presented in Nature  547, 28–331 (2017) and Donnelly et al., NJP 
(2018). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The 3D magnetisation vector field is reconstructed using the script 
    mag_tomo_recons.m
During the reconstruction, the magnetic structure is constrained to the 
magnetic material using a mask. This mask is calculated from the 
non-magnetic tomogram of the sample. Either the filename can be provided in 
the variable:
    nonmag_tomo_file
or, if this variable is left empty, the non-magnetic tomogram is calculated 
as a first step in the reconstruction, calling the function:
    non_mag_tomo.m

Once the mask is defined, the magnetic reconstruction is performed using a 
gradient-based optimisation procedure as described in Donnelly et al., NJP 
(2018). The gradients and error metric are defined in 
    grad_3D_recons_mag_mask.m

The following parameters should be filled in:
%%
Max number of iterations for the magnetic reconstruction - 
the more, the better. ~5-10 for a quick check, ~50+ for a more correct 
reconstruction (check convergence of error metric):
    it_max  
%%
Max number of iterations for the non-magnetic reconstruction (if no 
tomogram provided):
    it_max_mask 
%%
Choose whether you want to apply a mask to constrain the location of the 
magnetic material? 1=yes, 0=no. To avoid blurring of the edges, choose yes.
    applymask  
%%
Choose whether you want to save the reconstruction, and if so, where:
    save_recons
    save_file_name


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A quick guide to looking at the data in paraview:
    1. Import the data
    2. Threshold using the density tomogram to only select the magnetic material
To plot domain walls (e.g. mz=0):
    3. Use calculator to define mz.
    4. Use the contour function to define mz=0.
To plot vortex cores (e.g. in-plane vortices):
    5. Use calculator to define (mx_2+my_2)
    6. Use the contour function to define mx_2+my_2=(small value).
To plot a magnetic structure:
    7. Use a clip to select a part of the structure, e.g. a sphere or a box
       for a slice.
    8. Use the glyph function to represent the magnetisation by arrows.
OR
    9. Use the streamlines function to represent lines of constant 
       magnetisation by streamlines. The tube function can be used 
       afterwards for publication-quality figures
To export a high resolution figure:
    10. File -> Save Screenshot -> Choose number of pixels -> save.
