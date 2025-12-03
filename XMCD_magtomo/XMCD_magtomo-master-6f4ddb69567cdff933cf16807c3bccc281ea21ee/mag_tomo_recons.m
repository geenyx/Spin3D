%% Magnetic tomography with iterative optimisation

% Any work making use of these codes should refer to, and cite, the 
% following paper:
% "Tomographic reconstruction of a three-dimensional magnetization 
% vector field", Donnelly et al., NJP 2018.

% From a dataset measured with one set of circularly polarised light, 
% we recover three components of the magnetisation mx, my and mz
% We use arbitrary projections and simulataneously reconstruct all three
% components of the magnetisation

clear
addpath utilities/

% For the worked example provided, projections can be calculated of 
% micromagnetic simulations in mag_tomo_calculate_projections.m

% Here we define the base path. For the simulations, can choose to simulate
% a smaller cropped structure to save time. Choose cropped or full 
% projections here:
cropping = 1; % 1=cropped projections, 0=full

if cropping
    basepath = '/das/work/p16/p16766/example_code_run/Cropped_mag_structure/';
else
    basepath = '/das/work/p16/p16766/example_code_run/Full_mag_structure/';
end
%% Load the projections

% As calculated, define the projections as follows:

num_tilt_angles = 2; % how many tilt axes are used
theta_num = 180; % how many angles do you measure around 360 degrees?


load(sprintf('%smag_projectionsABS360_crop_%daxes_%dtheta.mat',basepath,num_tilt_angles,theta_num));


%% Define the parameters for magnetic tomography:


% provide the non-magnetic reconstruction here:
nonmag_tomo_file = '/das/work/p16/p16766/example_code_run/Cropped_mag_structure/non_mag_tomo_it_10.mat'; %[];
% if not available, leave the filename empty and it will be calculated 
% from the single polarisation projections.

it_max = 10; % Max number of iterations for the magnetic reconstruction
it_max_mask = 10; % Max number of iterations for the non-magnetic reconstruction (if no tomo provided)
applymask = 1; % Do you want to apply a mask to constrain the location of the magnetic material? 1=yes, 0=no.

save_recons = 1; % Save the reconstruction


theta_max = ceil(thetas(end));

% Define the filename for saving the reconstructions:
save_file_name = sprintf('%smag_recon_angles_1to%d_num_tilt_axes_%d_thetas_%d_it_%d.mat',basepath,theta_max,num_tilt_angles,theta_num,it_max);

%% First non-magnetic tomography:
% the non-magnetic tomogram will be used to define a mask for the magnetic
% reconstruction. If it's not available, it will be reconstructed.


if isempty(nonmag_tomo_file)
    
    [abs_out] = non_mag_tomo(proj_out_mag,R,it_max_mask); % this could probably be replaced by filter back projection...
    
    % save this for future use:
    if cropping
        save(sprintf('%sCropped_mag_structure/non_mag_tomo_it_%d.mat',basepath,it_max_mask),'abs_out')
    else
        save(sprintf('%sFull_mag_structure/non_mag_tomo_it_%d.mat',basepath,it_max_mask),'abs_out')
    end
    
else
   
    load(nonmag_tomo_file); % the tomogram should be called abs_out
    
end

% Define the mask of the material:
% For the micromagnetic simulations provided, the non-magnetic signal = 1
% Therefore the mask is defined by:
mask = abs_out>0.5; % adjust value depending on the absorption of your material.


% check whether it makes sense:

figure(1)
subplot(1,2,1)
imagesc(abs_out(:,:,ceil(size(abs_out,3)/2)))
axis xy equal tight
title('Non-mag tomogram')
subplot(1,2,2)
imagesc(mask(:,:,ceil(size(abs_out,3)/2)))
axis xy equal tight
title('Mask')



%% Magnetic tomography

tic

[mx_out, my_out, mz_out,abs_out] = mag_tomo(proj_out_mag,R,it_max,abs_out,mask,applymask);
toc


% Save the reconstruction
if save_recons
    save(save_file_name,'mx_out','my_out','mz_out','abs_out','mask');
end



%% Some figures of the reconstruction:

ellipse_colour = [-1.5 1.5];
slice = 50;
figure(1002);
subplot(4,1,1)
imagesc(squeeze(mx_out(:,:,slice)));
title('mx '); 
axis equal tight xy;
colorbar;
 caxis(ellipse_colour)
colormap jet
subplot(4,1,2)
imagesc(squeeze(my_out(:,:,slice)));
title('my '); 
axis equal tight xy;
colorbar;
 caxis(ellipse_colour)
colormap jet
subplot(4,1,3)
imagesc(squeeze(mz_out(:,:,slice)));
title('mz '); 
axis equal tight xy;
colorbar;
 caxis(ellipse_colour)
colormap jet
subplot(4,1,4)
imagesc(squeeze(abs_out(:,:,slice)));
title('abs'); 
axis equal tight xy;
colorbar;
 caxis(ellipse_colour)
colormap jet

%% Write a vtk file, compatible with paraview for viewing data


X = [1:size(mx_out,1)];
Y = [1:size(mx_out,2)];
Z = [1:size(mx_out,3)];
[x,y,z] = meshgrid(X,Y,Z);


vtkwrite(sprintf('%smag_recons_%daxes_%dtheta.vtk',basepath,num_tilt_angles,theta_num),'structured_grid',x,y,z,'vectors','magnetisation',mx_out,my_out,mz_out,'scalars','density',abs_out);%,'scalars','tomo_recons',tomo_recons);



% Copyright 2018, ETH Zurich, Claire Donnelly, Manuel Guizar-Sicairos

% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:

% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
