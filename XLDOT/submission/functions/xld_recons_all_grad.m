function [E,grad] = xld_recons_all_grad(opt_inputs,params,proj_xmld,pol_vec)
    R = params.R;
    % p = params.p;

    cfg=params.cfg;
    vectors=params.vectors;
    split=params.split;

    mask = params.mask;

    mx = opt_inputs(:,:,:,1);
    my = opt_inputs(:,:,:,2);
    mz = opt_inputs(:,:,:,3);

    clear abs_vect mx_vect my_vect mz_vect

    %% First calculate the projections of the guess  structure:
    proj_xx=astra.Ax_partial(mx.*mx, cfg, vectors, split);
    proj_xy=astra.Ax_partial(mx.*my, cfg, vectors, split);
    proj_xz=astra.Ax_partial(mx.*mz, cfg, vectors, split);
    proj_yy=astra.Ax_partial(my.*my, cfg, vectors, split);
    proj_yz=astra.Ax_partial(my.*mz, cfg, vectors, split);
    proj_zz=astra.Ax_partial(mz.*mz, cfg, vectors, split);
    proj_guess_xmld=zeros(size(proj_xx,1),size(proj_xx,2),size(proj_xx,3));
    
    for ii = 1:cfg.iProjAngles
        RR = R(:,:,ii);
        pol = pol_vec(ii);
        proj_guess_xmld(:,:,ii) = ((proj_xx(:,:,ii)  *(RR(1,1)*cosd(pol)+RR(2,1)*sind(pol))^2 + ...
                                  proj_yy(:,:,ii)  *(RR(1,2)*cosd(pol)+RR(2,2)*sind(pol))^2 + ...
                                  proj_zz(:,:,ii)  *(RR(1,3)*cosd(pol)+RR(2,3)*sind(pol))^2 + ...
                                  2*proj_xy(:,:,ii)*(RR(1,1)*cosd(pol)+RR(2,1)*sind(pol))*(RR(1,2)*cosd(pol)+RR(2,2)*sind(pol)) + ...
                                  2*proj_xz(:,:,ii)*(RR(1,1)*cosd(pol)+RR(2,1)*sind(pol))*(RR(1,3)*cosd(pol)+RR(2,3)*sind(pol)) + ...
                                  2*proj_yz(:,:,ii)*(RR(1,2)*cosd(pol)+RR(2,2)*sind(pol))*(RR(1,3)*cosd(pol)+RR(2,3)*sind(pol))));
    end

    %% A check to make sure the arrays match:
    if ~isequal(size(proj_xmld),size(proj_guess_xmld))
        disp('Oops, the projection and the guess are different sizes - please check');
        keyboard
    end
    %% Calculate error metric:

    proj_diff=single(proj_guess_xmld - proj_xmld);
    E = sum(proj_diff(:).^2);

    %% Calculate the gradient of the error metric:


    proj_diff_xx=single(zeros(size(proj_diff,1),size(proj_diff,2),size(proj_diff,3)));
    proj_diff_xy=single(zeros(size(proj_diff_xx)));
    proj_diff_xz=single(zeros(size(proj_diff_xx)));
    proj_diff_yx=single(zeros(size(proj_diff_xx)));
    proj_diff_yy=single(zeros(size(proj_diff_xx)));
    proj_diff_yz=single(zeros(size(proj_diff_xx)));
    proj_diff_zx=single(zeros(size(proj_diff_xx)));
    proj_diff_zy=single(zeros(size(proj_diff_xx)));
    proj_diff_zz=single(zeros(size(proj_diff_xx)));

    for ii = 1:cfg.iProjAngles
        RR = R(:,:,ii);       
        cp = cosd(pol_vec(ii));
        sp = sind(pol_vec(ii));

        proj_diff_xx(:,:,ii)=proj_diff(:,:,ii) * (RR(1,1)*cp + RR(2,1)*sp) * ...
            (RR(1,1)*cp + RR(2,1)*sp);
        proj_diff_xy(:,:,ii)=proj_diff(:,:,ii) * (RR(1,1)*cp + RR(2,1)*sp) * ...
            (RR(1,2)*cp + RR(2,2)*sp); 
        proj_diff_xz(:,:,ii)=proj_diff(:,:,ii) * (RR(1,1)*cp + RR(2,1)*sp) * ...
            (RR(1,3)*cp + RR(2,3)*sp);

        proj_diff_yx(:,:,ii)=proj_diff(:,:,ii) * (RR(1,2)*cp + RR(2,2)*sp) * ...
            (RR(1,1)*cp + RR(2,1)*sp);
        proj_diff_yy(:,:,ii)=proj_diff(:,:,ii) * (RR(1,2)*cp + RR(2,2)*sp) * ...
            (RR(1,2)*cp + RR(2,2)*sp);
        proj_diff_yz(:,:,ii)=proj_diff(:,:,ii) * (RR(1,2)*cp + RR(2,2)*sp) * ...
            (RR(1,3)*cp + RR(2,3)*sp);    

        proj_diff_zx(:,:,ii)=proj_diff(:,:,ii) * (RR(1,3)*cp + RR(2,3)*sp) * ...
            (RR(1,1)*cp + RR(2,1)*sp);
        proj_diff_zy(:,:,ii)=proj_diff(:,:,ii) * (RR(1,3)*cp + RR(2,3)*sp) * ...
            (RR(1,2)*cp + RR(2,2)*sp);
        proj_diff_zz(:,:,ii)=proj_diff(:,:,ii) * (RR(1,3)*cp + RR(2,3)*sp) * ...
            (RR(1,3)*cp + RR(2,3)*sp);           

    end

    proj_diff_xx = single(proj_diff_xx);
    proj_diff_xy = single(proj_diff_xy);
    proj_diff_xz = single(proj_diff_xz);

    proj_diff_yx = single(proj_diff_yx);
    proj_diff_yy = single(proj_diff_yy);
    proj_diff_yz = single(proj_diff_yz);

    proj_diff_zx = single(proj_diff_zx);
    proj_diff_zy = single(proj_diff_zy);
    proj_diff_zz = single(proj_diff_zz);

    grad_xx = 2*astra.Atx_partial(proj_diff_xx, cfg, vectors, split);
    grad_xy = 2*astra.Atx_partial(proj_diff_xy, cfg, vectors, split);
    grad_xz = 2*astra.Atx_partial(proj_diff_xz, cfg, vectors, split);
    grad_yx = 2*astra.Atx_partial(proj_diff_yx, cfg, vectors, split);
    grad_yy = 2*astra.Atx_partial(proj_diff_yy, cfg, vectors, split);
    grad_yz = 2*astra.Atx_partial(proj_diff_yz, cfg, vectors, split);
    grad_zx = 2*astra.Atx_partial(proj_diff_zx, cfg, vectors, split);
    grad_zy = 2*astra.Atx_partial(proj_diff_zy, cfg, vectors, split);
    grad_zz = 2*astra.Atx_partial(proj_diff_zz, cfg, vectors, split);

    grad_x  = (grad_xx.*mx + grad_xy.*my + grad_xz.*mz);
    grad_y  = (grad_yx.*mx + grad_yy.*my + grad_yz.*mz);
    grad_z  = (grad_zx.*mx + grad_zy.*my + grad_zz.*mz);

    grad_x = grad_x.*mask;
    grad_y = grad_y.*mask;
    grad_z = grad_z.*mask;

    grad = cat(4, grad_x, grad_y, grad_z);

end