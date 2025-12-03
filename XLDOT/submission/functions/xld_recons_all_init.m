function [mx_out, my_out, mz_out, E] = xld_recons_all_init(xmld_proj,R,it_max,mask,params,pol_vec,mx_init,my_init,mz_init)
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

    opt_inputs = cat(4,mx,my,mz);     % A vector containing all values of the input to the reconstruction, both mx, my and mz

    tic
        %%%%% Optimization with gradient %%%%%%
        itmax=[it_max]; %   Maximum number of iterations (empty for default = 50)
        ftol=[1e-7];    %   Relative function tolerance (empty for default = 1e-3)
        xtol=[1e-7];    %   Absolute solution tolerance (empty for default = 1e-3)
        [opt_out, E] = cgmin1('xld_recons_all_grad',opt_inputs,itmax,ftol,xtol,params,xmld_proj,pol_vec);
    toc
        mx_out = opt_out(:,:,:,1);
        my_out = opt_out(:,:,:,2);
        mz_out = opt_out(:,:,:,3);
end