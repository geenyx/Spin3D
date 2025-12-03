% function [E,grad] = grad_3D_recons(opt_inputs,params,proj_mag_out)
%
% Calculates the error metric and gradient for a non-magnetic tomogram 
% measured with single circular polarisation XMCD projections.
%
% opt_inputs - contains the input "guess" of the variables that will be
% reconstructed. Each variable is a column of opt_inputs - mx, my, mz, and
% abs, the non-magnetic tomo.
%
% params - contains useful parameters for the magnetic reconstruction that
% are defined in mag_tomo.m. These include the size of the tomogram, the
% size of the set of rotation matrices, the mask and "applymask".
%
% proj_mag_out - a stack of the measured projections


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

function [E,grad] = grad_3D_recons(opt_inputs,params,proj_mag_out)

mynargout = nargout;
numx = params.numx;                                 % Size of tomo in x
numy = params.numy;                                 % Size of tomo in y
numz = params.numz;                                 % Size of tomo in z
Rx = params.Rx;                                 % Size of tomo in x
Ry = params.Ry;                                 % Size of tomo in y
Rz = params.Rz;                                 % Size of tomo in z
R = params.R;
p = params.p;

abs_vect = opt_inputs(:,1);
ABS = reshape(abs_vect,numx,numy,numz);


%% First calculate the projections of the guess  structure:


x = [1:size(ABS,1)]-ceil(size(ABS,1)/2);
y = [1:size(ABS,3)]-ceil(size(ABS,3)/2);
z = [1:size(ABS,2)]-ceil(size(ABS,2)/2);

[X,Z,Y] = meshgrid(x,z,y);

for ii = 1:size(R,3)
[proj_out_guess(:,:,ii),xout,yout] = arb_projection(ABS, X,Y,Z, squeeze(R(:,:,ii)), p);%, xout, yout );
end

%% Calculate error metric:

proj_diff = (proj_out_guess-proj_mag_out); %difference between the estimated and measured projections
proj_diff_sq = (proj_out_guess-proj_mag_out).^2; %take the square

E = sum(proj_diff_sq(:)); % calculate error metric

%% Calculate the gradient of the error metric:

grad_recons = zeros(size(X));


if mynargout > 1    % Only calculate the gradient if the function asks for more than one output
    for ii = 1:size(R,3)
          back_proj = arb_back_projection( proj_diff(:,:,ii), xout, yout, X, Y, Z, (R(:,:,ii)), p)   ;
          grad_recons = grad_recons+2*back_proj;
    end
    


grad = [grad_recons(:)];
end
end
