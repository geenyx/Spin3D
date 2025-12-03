function [mx_out, my_out, mz_out,E] = xld_recons_diff_init_state(xld_proj,R,it_max,mask,params,pol_H,pol_V,pol_purity,mx_init,my_init,mz_init)
%% Define the parameters for tomography

% define the parameters of the projections:
p.volume_upsampling = 0;
p.filter_2D = 3;           % Lower values give you compromise between resolution and artifacts
p.method = 'bilinear';  

cfg = params.cfg;

mx = single(mx_init.*mask);
my = single(my_init.*mask);
mz = single(mz_init.*mask);

%%
% Define parameters needed for reconstruction
params.R = R;
params.p = p;                           % Parameters for projections

%% Now the magnetic structure, with a mask

% applymask = 1;

    params.mask=mask;             % apply a mask with the electron density

    opt_inputs = cat(4,mx,my,mz);    % A vector containing all values of the input to the reconstruction, both mx, my and mz   

tic
    %%%%% Optimization with gradient %%%%%%
    itmax=it_max; %   Maximum number of iterations (empty for default = 50)
    ftol=1e-7;    %   Relative function tolerance (empty for default = 1e-3)
    xtol=1e-7;    %   Absolute solution tolerance (empty for default = 1e-3)
    [opt_out, E] = cgmin1('xld_recons_diff_grad',opt_inputs,itmax,ftol,xtol,params,xld_proj,pol_H,pol_V,pol_purity); %
toc
    mx_out = opt_out(:,:,:,1);
    my_out = opt_out(:,:,:,2);
    mz_out = opt_out(:,:,:,3);
%     abs_out = reshape(opt_out1(:,4),size(mz,1),size(mz,2),size(mz,3));
    
    %save('/das/work/p16/p16766/3D_recons/Cropped_mag_structure/mag_recon_1_180_0_30_-30_it_500_mask.mat','mx_out','my_out','mz_out');