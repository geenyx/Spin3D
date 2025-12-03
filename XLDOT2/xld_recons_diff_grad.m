function [E,grad] = xld_recons_diff_grad(opt_inputs,params,proj_xld,pol_H,pol_V,pol_purity)
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
guess_H=zeros(size(proj_xx,1),size(proj_xx,2),size(proj_xx,3));
guess_V=zeros(size(proj_xx,1),size(proj_xx,2),size(proj_xx,3));

for ii = 1:cfg.iProjAngles
    RR = R(:,:,ii);
    guess_H(:,:,ii) = proj_xx(:,:,ii)  *(RR(1,1)*cosd(pol_H)+RR(2,1)*sind(pol_H))^2 + ...
                      proj_yy(:,:,ii)  *(RR(1,2)*cosd(pol_H)+RR(2,2)*sind(pol_H))^2 + ...
                      proj_zz(:,:,ii)  *(RR(1,3)*cosd(pol_H)+RR(2,3)*sind(pol_H))^2 + ...
                      2*proj_xy(:,:,ii)*(RR(1,1)*cosd(pol_H)+RR(2,1)*sind(pol_H))*(RR(1,2)*cosd(pol_H)+RR(2,2)*sind(pol_H)) + ...
                      2*proj_xz(:,:,ii)*(RR(1,1)*cosd(pol_H)+RR(2,1)*sind(pol_H))*(RR(1,3)*cosd(pol_H)+RR(2,3)*sind(pol_H)) + ...
                      2*proj_yz(:,:,ii)*(RR(1,2)*cosd(pol_H)+RR(2,2)*sind(pol_H))*(RR(1,3)*cosd(pol_H)+RR(2,3)*sind(pol_H));
    guess_V(:,:,ii) = (proj_xx(:,:,ii) *(RR(1,1)*cosd(pol_V)+RR(2,1)*sind(pol_V))^2 + ...
                      proj_yy(:,:,ii)  *(RR(1,2)*cosd(pol_V)+RR(2,2)*sind(pol_V))^2 + ...
                      proj_zz(:,:,ii)  *(RR(1,3)*cosd(pol_V)+RR(2,3)*sind(pol_V))^2 + ...
                      2*proj_xy(:,:,ii)*(RR(1,1)*cosd(pol_V)+RR(2,1)*sind(pol_V))*(RR(1,2)*cosd(pol_V)+RR(2,2)*sind(pol_V)) + ...
                      2*proj_xz(:,:,ii)*(RR(1,1)*cosd(pol_V)+RR(2,1)*sind(pol_V))*(RR(1,3)*cosd(pol_V)+RR(2,3)*sind(pol_V)) + ...
                      2*proj_yz(:,:,ii)*(RR(1,2)*cosd(pol_V)+RR(2,2)*sind(pol_V))*(RR(1,3)*cosd(pol_V)+RR(2,3)*sind(pol_V)))./pol_purity;
end
guess_xld = guess_H - guess_V;
clear proj_xx proj_xy proj_xz proj_yy proj_yz proj_zz guess_H guess_V
%% A check to make sure the arrays match:
if ~isequal(size(proj_xld),size(guess_xld))
    disp('Oops, the projection and the guess are different sizes - please check');
    keyboard
end

%% Calculate error metric:
proj_diff=single(guess_xld - proj_xld);
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
    
    cp_H = cosd(pol_H);
    sp_H = sind(pol_H);
    
    cp_V = cosd(pol_V);
    sp_V = sind(pol_V);
    
    for ii = 1:cfg.iProjAngles
        RR = R(:,:,ii);       
        proj_diff_xx(:,:,ii)=proj_diff(:,:,ii) * ((RR(1,1)*cp_H + RR(2,1)*sp_H) * ...
            (RR(1,1)*cp_H + RR(2,1)*sp_H) - (RR(1,1)*cp_V + RR(2,1)*sp_V) * ...
            (RR(1,1)*cp_V + RR(2,1)*sp_V));
        
        proj_diff_xy(:,:,ii)=proj_diff(:,:,ii) * ((RR(1,1)*cp_H + RR(2,1)*sp_H) * ...
            (RR(1,2)*cp_H + RR(2,2)*sp_H) - (RR(1,1)*cp_V + RR(2,1)*sp_V) * ...
            (RR(1,2)*cp_V + RR(2,2)*sp_V)); 
        
        proj_diff_xz(:,:,ii)=proj_diff(:,:,ii) * ((RR(1,1)*cp_H + RR(2,1)*sp_H) * ...
            (RR(1,3)*cp_H + RR(2,3)*sp_H) - (RR(1,1)*cp_V + RR(2,1)*sp_V) * ...
            (RR(1,3)*cp_V + RR(2,3)*sp_V));
        
        
        proj_diff_yx(:,:,ii)=proj_diff(:,:,ii) * ((RR(1,2)*cp_H + RR(2,2)*sp_H) * ...
            (RR(1,1)*cp_H + RR(2,1)*sp_H) - (RR(1,2)*cp_V + RR(2,2)*sp_V) * ...
            (RR(1,1)*cp_V + RR(2,1)*sp_V));
        proj_diff_yy(:,:,ii)=proj_diff(:,:,ii) * ((RR(1,2)*cp_H + RR(2,2)*sp_H) * ...
            (RR(1,2)*cp_H + RR(2,2)*sp_H) - (RR(1,2)*cp_V + RR(2,2)*sp_V) * ...
            (RR(1,2)*cp_V + RR(2,2)*sp_V));
        proj_diff_yz(:,:,ii)=proj_diff(:,:,ii) * ((RR(1,2)*cp_H + RR(2,2)*sp_H) * ...
            (RR(1,3)*cp_H + RR(2,3)*sp_H) - (RR(1,2)*cp_V + RR(2,2)*sp_V) * ...
            (RR(1,3)*cp_V + RR(2,3)*sp_V));    
        
        proj_diff_zx(:,:,ii)=proj_diff(:,:,ii) * ((RR(1,3)*cp_H + RR(2,3)*sp_H) * ...
            (RR(1,1)*cp_H + RR(2,1)*sp_H) - (RR(1,3)*cp_V + RR(2,3)*sp_V) * ...
            (RR(1,1)*cp_V + RR(2,1)*sp_V));
        proj_diff_zy(:,:,ii)=proj_diff(:,:,ii) * ((RR(1,3)*cp_H + RR(2,3)*sp_H) * ...
            (RR(1,2)*cp_H + RR(2,2)*sp_H) - (RR(1,3)*cp_V + RR(2,3)*sp_V) * ...
            (RR(1,2)*cp_V + RR(2,2)*sp_V));
        proj_diff_zz(:,:,ii)=proj_diff(:,:,ii) * ((RR(1,3)*cp_H + RR(2,3)*sp_H) * ...
            (RR(1,3)*cp_H + RR(2,3)*sp_H) - (RR(1,3)*cp_V + RR(2,3)*sp_V) * ...
            (RR(1,3)*cp_V + RR(2,3)*sp_V));           
    end
    
    grad_xx = 2*astra.Atx_partial(proj_diff_xx, cfg, vectors, split);
    grad_xy = 2*astra.Atx_partial(proj_diff_xy, cfg, vectors, split);
    grad_xz = 2*astra.Atx_partial(proj_diff_xz, cfg, vectors, split);
    grad_yx = 2*astra.Atx_partial(proj_diff_yx, cfg, vectors, split);
    grad_yy = 2*astra.Atx_partial(proj_diff_yy, cfg, vectors, split);
    grad_yz = 2*astra.Atx_partial(proj_diff_yz, cfg, vectors, split);
    grad_zx = 2*astra.Atx_partial(proj_diff_zx, cfg, vectors, split);
    grad_zy = 2*astra.Atx_partial(proj_diff_zy, cfg, vectors, split);
    grad_zz = 2*astra.Atx_partial(proj_diff_zz, cfg, vectors, split);

    clear proj_diff_xx proj_diff_xy proj_diff_xz proj_diff_yx proj_diff_yy
    clear proj_diff_yz proj_diff_zx proj_diff_zy proj_diff_zz
    
    grad_x  = grad_xx.*mx + grad_xy.*my + grad_xz.*mz;
    grad_y  = grad_yx.*mx + grad_yy.*my + grad_yz.*mz;
    grad_z  = grad_zx.*mx + grad_zy.*my + grad_zz.*mz;
    
    grad_x = grad_x.*mask;
    grad_y = grad_y.*mask;
    grad_z = grad_z.*mask;

    grad = cat(4, grad_x, grad_y, grad_z);
end
