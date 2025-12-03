% function [mx_out, my_out, mz_out,abs_out] = mag_tomo(proj_out_mag,R,it_max,abs_out,mask,applymask)
% 
% Returns the three components of the magnetisation and the non-magnetic
% tomogram for an input of single circular polarisation projections.
% 
% proj_out_mag - a stack of the single circular polarisation projections.
% 
% R - a stack of the rotation matrices that define the orientation of the
% sample for each projection.
% 
% it_max - the maximum number of iterations for the reconstruction of 
% the magnetic tomogram.
% 
% abs_out - the non-magnetic tomogram.
% 
% mask - a logical mask for the location of the magnetic material that is
% calculated from the non-magnetic tomogram. 1 = material, 0=air.
% 
% applymask - choose whether the magnetisation is constrained to the mask
% during the reconstruction.


% Copyright 2018, ETH Zurich, Claire Donnelly, Manuel Guizar-Sicairos
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
% 
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


% Any work making use of these codes should refer to, and cite, the 
% following paper:
% "Tomographic reconstruction of a three-dimensional magnetization 
% vector field", Donnelly et al., NJP 2018.


function [mx_out, my_out, mz_out,abs_out] = mag_tomo(proj_out_mag,R,it_max,abs_out,mask,applymask)


%% Define the parameters for tomography

% define the parameters of the projections:
p.volume_upsampling = 0;
p.filter_2D = 3;           % Lower values give you compromise between resolution and artifacts
p.method = 'bilinear';   


% initial guess is zeros:

mx = zeros(size(proj_out_mag,2),size(proj_out_mag,2),size(proj_out_mag,1));
my = zeros(size(proj_out_mag,2),size(proj_out_mag,2),size(proj_out_mag,1)); 
mz = zeros(size(proj_out_mag,2),size(proj_out_mag,2),size(proj_out_mag,1));
ABS = zeros(size(proj_out_mag,2),size(proj_out_mag,2),size(proj_out_mag,1));

%% Define parameters needed for reconstruction
params.numx = size(mz,1);               % Size of mz_x
params.numy = size(mz,2);               % Size of mz_y
params.numz = size(mz,3);               % Size of mz_z
params.Rx = size(R,1);                  % Size of R_x
params.Ry = size(R,2);                  % Size of R_y
params.Rz = size(R,3);                  % Size of R_z
params.R = R;
params.p = p;                           % Parameters for projections

%% Now the magnetic structure, with a mask

params.applymask=applymask;             % apply a mask with the electron density
params.mask=mask;             % apply a mask with the electron density

opt_inputs = [mx(:) my(:) mz(:) abs_out(:)];    % A vector containing all values of the input to the reconstruction, both mx, my and mz    


%%%%% Optimization with gradient %%%%%%
itmax=[it_max]; %   Maximum number of iterations (empty for default = 50)
ftol=[1e-5];    %   Relative function tolerance (empty for default = 1e-3)
xtol=[1e-5];    %   Absolute solution tolerance (empty for default = 1e-3)
opt_out1 = cgmin1('grad_3D_recons_mag_mask',opt_inputs,itmax,ftol,xtol,params, proj_out_mag); %


mx_out = reshape(opt_out1(:,1),size(mz,1),size(mz,2),size(mz,3));
my_out = reshape(opt_out1(:,2),size(mz,1),size(mz,2),size(mz,3));
mz_out = reshape(opt_out1(:,3),size(mz,1),size(mz,2),size(mz,3));
abs_out = reshape(opt_out1(:,4),size(mz,1),size(mz,2),size(mz,3));

end