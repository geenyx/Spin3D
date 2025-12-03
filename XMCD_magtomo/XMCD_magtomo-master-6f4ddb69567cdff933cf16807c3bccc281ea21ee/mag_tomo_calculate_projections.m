%% Magnetic tomography with iterative optimisation


% Any work making use of these codes should refer to, and cite, the 
% following paper:
% "Tomographic reconstruction of a three-dimensional magnetization 
% vector field", Donnelly et al., NJP 2018.


% This script is used to produce single projection X-ray projections, with 
% XMCD magnetic signal, from a micromagnetic simulation.


clear
addpath utilities

basepath = '/das/work/p16/p16766/example_code_run/';


data = load(sprintf('%scylinder_end.mat',basepath));

cropping=1; % take a smaller section of the pillar to reduce calculation times.

if cropping
    
    % crop the structure for a quicker calculation for testing
    MX_=data.MX(:,:,31:80);
    MY_=data.MY(:,:,31:80);
    MZ_=data.MZ(:,:,31:80);
    ABS_=data.ABS(:,:,31:80);
    pixsize=data.pixsize;

    %%% pad out the structure with zeros so that the projections are not cropped
    MX = zeros(size(MX_,1)+30,size(MX_,2)+30,size(MX_,3)+50);
    MY = zeros(size(MX_,1)+30,size(MX_,2)+30,size(MX_,3)+50);
    MZ = zeros(size(MX_,1)+30,size(MX_,2)+30,size(MX_,3)+50);
    ABS = zeros(size(MX_,1)+30,size(MX_,2)+30,size(MX_,3)+50);

    MX(16:end-15,16:end-15,26:end-25)=MX_(:,:,:);
    MY(16:end-15,16:end-15,26:end-25)=MY_(:,:,:);
    MZ(16:end-15,16:end-15,26:end-25)=MZ_(:,:,:);
    ABS(16:end-15,16:end-15,26:end-25)=ABS_(:,:,:);
    
else 
    
    MX_=data.MX(:,:,:);
    MY_=data.MY(:,:,:);
    MZ_=data.MZ(:,:,:);
    ABS_=data.ABS(:,:,:);
    pixsize=data.pixsize;

    %%% pad out the structure with zeros so that the projections are not cropped
    MX = zeros(size(MX_,1)+90,size(MX_,2)+90,size(MX_,3)+40);
    MY = zeros(size(MX_,1)+90,size(MX_,2)+90,size(MX_,3)+40);
    MZ = zeros(size(MX_,1)+90,size(MX_,2)+90,size(MX_,3)+40);
    ABS = zeros(size(MX_,1)+90,size(MX_,2)+90,size(MX_,3)+40);

    MX(46:end-45,46:end-45,21:end-20)=MX_(:,:,:);
    MY(46:end-45,46:end-45,21:end-20)=MY_(:,:,:);
    MZ(46:end-45,46:end-45,21:end-20)=MZ_(:,:,:);
    ABS(46:end-45,46:end-45,21:end-20)=ABS_(:,:,:);

end



%% Define the geometry
% Here the geometry is based on the dual axis tomography demonstrated in
% Nature 547, 328–331 (2017), and in NJP ???

% tilt angles: (in degrees)
thetas = [0:2:360-1e-3]; % Angles measured around rotation axis (tomographic measurement)
tilt_angles = [0 30]; % Tilt angles for different tomo measurements

num_tilt_angles=size(tilt_angles,2);
theta_num=size(thetas,2);


%% First define the arbitrary projections:

% Each projection will be defined by a rotation matrix, which will be the
% product of two rotation matrices in this case of dual axis tomography.
% The structure is first rotated around the z (2nd index), Rot_z, and then 
% rotating around y, Rot_y. 

R = zeros(3,3,size(tilt_angles,1)*size(thetas,1));
ii=0;
for tilt_angle = tilt_angles
for theta =[thetas]
    ii = ii+1;
    % tilt about the z axis
    Rot_z=[cosd(tilt_angle) sind(tilt_angle) 0
            -sind(tilt_angle)  cosd(tilt_angle) 0
            0 0 1];
    % measure a tomogram about the y rotation axis
    Rot_y = [cosd(theta) 0 sind(theta) 
            0 1 0
            -sind(theta) 0 cosd(theta)];
    R(:,:,ii) = Rot_y*Rot_z;
end
end

%% Calculate the projections:


x = [1:size(MX,1)]-ceil(size(MX,1)/2);
y = [1:size(MX,3)]-ceil(size(MX,3)/2);
z = [1:size(MX,2)]-ceil(size(MX,2)/2);

[X,Z,Y] = meshgrid(x,z,y);

% define the parameters of the projection:
p.volume_upsampling = 0;
p.filter_2D = 3;           % Lower values give you compromise between resolution and artifacts
p.method = 'bilinear';   

ii=0;
for tilt_angle = tilt_angles
for theta =[thetas]
ii = 1+ii;
R_ = (squeeze(R(:,:,ii)));
MY_rotate = (R_(3,1)*MX + R_(3,2)*MZ +R_(3,3)*MY)+ABS;

[ proj_out_mag(:,:,ii), xout, yout ] = arb_projection(MY_rotate, X,Y,Z, squeeze(R(:,:,ii)), p);

end
end


%% Look at some figures:

% Choose ii to be for angles 0:180, then the magnetic "XMCD "projection 
% will be calculated from the two antiparallel projections.

ii=1;
figure(11)
subplot(1,2,1)
imagesc(proj_out_mag(:,:,ii))
axis xy equal tight
colorbar
title('Magnetic projection')
subplot(1,2,2)
imagesc(proj_out_mag(:,:,ii)-flip(proj_out_mag(:,:,ii+theta_num/2),2))
axis xy equal tight
title('Magnetic component')
colorbar


%% sanity check to check magnetic structure
% for the simulations, you can have a look and make sure your projections
% make sense.

slice = 30;
figure(2)
subplot(3,1,1)
imagesc(MX(:,:,slice))
axis xy equal tight
caxis([-1.5 1.5])
colorbar
title('mx')
subplot(3,1,2)
imagesc(MY(:,:,slice))
axis xy equal tight
caxis([-1.5 1.5])
colorbar
title('my')
subplot(3,1,3)
imagesc(MZ(:,:,slice))
axis xy equal tight
caxis([-1.5 1.5])
colorbar
title('mz')
%  caxis([-30 30])
colormap(jet)

%% Save the projections

if cropping
    save(sprintf('%sCropped_mag_structure/mag_projectionsABS360_crop_%daxes_%dtheta.mat',basepath,num_tilt_angles,theta_num),'thetas','tilt_angles','proj_out_mag','R');
else
    save(sprintf('%sFull_mag_structure/mag_projections360_full_%daxes_%dtheta.mat',basepath,num_tilt_angles,theta_num),'thetas','tilt_angles','proj_out_mag','R');
end


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
