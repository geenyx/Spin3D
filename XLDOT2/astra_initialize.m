function [cfg, vectors] = astra_initialize(Npix,R)
% Creates configuration and vector variables that work with the ASTRA toolkit
% Angles are assumed in degrees and should be sorted to maximise performance
% (better use of texture memory)
assert(math.isint(Npix), 'Npix is not integer');

if isscalar(Npix)
    Npix(2) = Npix;
end
if length(Npix) == 2 && all(lamino_angle == 90)
    Npix(3) = size_projection(1);  % default behaviour is to have same number of layers in reconstruction and in laminography
elseif length(Npix) == 2 && any(lamino_angle ~= 90)
    error('All three dimensions of the volume size has to be specified for the laminograhy geometry')
end

%% Define cfg (configuration)
% Define the reconstruction volume in pixels using Npix
[cfg.iVolX,cfg.iVolY,cfg.iVolZ] = deal(Npix(1),Npix(2),Npix(3));
% Projection size (as seeon on detector, i.e. summing along index 1)
[cfg.iProjU,cfg.iProjV] = deal(Npix(2),Npix(3));
% Number of angles
cfg.iProjAngles = size(R,3);          
cfg.iRaysPerDet = 1;
cfg.iRaysPerDetDim = 1;
cfg.iRaysPerVoxelDim = 1;
cfg.pixel_scale = [1,1];           % px size on the detector [horizontal , vertical]
cfg.tilt_angle = 0;                % rotation in the detector plane
cfg.skewness_angle = 0;

%% Define vectors
% R is defined as in the arbitrary projection code by Manuel
% This integrates along z (2nd index) and so this code follows this
% convention.
R_conv=zeros(size(R));
convert_matrix =   [0 0 1
                    0 1 0
                    1 0 0];
for ii=1:size(R_conv,3)
    R_conv(:,:,ii)=R(:,:,ii)*convert_matrix;
end
% Flip x and z dimensions (index 1 and 3) for some reason...
%   maybe this is why the comment above says it inegrates along z?

% The astra documentation defines "vectors" as a stack of 4, 3d vectors:
% ray direction                                     (1,2,3)
% center of detector                                (4,5,6)
% vector from detector pixel (0,0) to (0,1)         (7,8,9)
% vector from detector pixel (0,0) to (1,0)         (10,11,12)

rotation_center = [Npix(3)/2, Npix(2)/2];  % around center of the projection
vectors = zeros(cfg.iProjAngles, 12);
for i = 1:cfg.iProjAngles
    vectors(i,1) = -R(3,1,i);
    vectors(i,2) = -R(3,3,i);
    vectors(i,3) = -R(3,2,i);

    vectors(i,4:6) = 0;

    vectors(i,7) = R(1,1,i)/cfg.pixel_scale(1);
    vectors(i,8) = R(1,3,i)/cfg.pixel_scale(1);
    vectors(i,9) = R(1,2,i)/cfg.pixel_scale(1);

    vectors(i,10) = R(2,1,i)/cfg.pixel_scale(2);
    vectors(i,11) = R(2,3,i)/cfg.pixel_scale(2);
    vectors(i,12) = R(2,2,i)/cfg.pixel_scale(2);

    vectors(i,4:6) = vectors(i,4:6) -(vectors(i,10:12)*(rotation_center(1))+vectors(i,7:9)*(rotation_center(2)));
end

%% Example from ASTRA
% % ray direction
% vectors(i,1) = sin(proj_geom.ProjectionAngles(i));
% vectors(i,2) = -cos(proj_geom.ProjectionAngles(i));
% vectors(i,3) = 0;
% 
% % center of detector
% vectors(i,4) = 0;
% vectors(i,5) = 0;
% vectors(i,6) = 0;
% 
% % vector from detector pixel (0,0) to (0,1)
% vectors(i,7) = cos(proj_geom.ProjectionAngles(i)) * proj_geom.DetectorSpacingX;
% vectors(i,8) = sin(proj_geom.ProjectionAngles(i)) * proj_geom.DetectorSpacingX;
% vectors(i,9) = 0;
% 
% % vector from detector pixel (0,0) to (1,0)
% vectors(i,10) = 0;
% vectors(i,11) = 0;
% vectors(i,12) = proj_geom.DetectorSpacingY;

% Focusing on the 3 x 3 (or rotation matrix equivalent vectors)
% looks like: [s -c  0] or its transpose [s  c  0]
%             [c  s  0]                  [-c s  0]
%             [0  0  1]                  [0  0  1]
% adding 90 deg to these:
% A :[c  s  0]    B: [c -s  0]
%    [-s c  0]       [s  c  0]
%    [0  0  1]       [0  0  1]
% which now look like the usual rotation matrices, rotating about z 
