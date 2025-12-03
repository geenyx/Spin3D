%% XLD reconstruction
clear
close all

addpath functions

%% Load dataset
load('reference_structures/smoothly_varying.mat')

%% Save path
path_name = 'output/smoothly_varying/dual_axis';

%%
pol = 0;

n_angles = 270;
tilts = [0,30];

R = eul2rot('zy',[0,30],linspace(0,179+1/3,n_angles));
Npix = size(mx);
pol_vec = [ones(1,n_angles)*pol, ones(1,n_angles)*pol,ones(1,n_angles)*(pol+90),ones(1,n_angles)*(pol+90)];

% Get projections
[xld_H,xld_V] = calculate_xld(mx,my,mz,R,pol,pol+90);

%% Reconstruct
R = cat(3,R,R);

[cfg, vectors] = astra_initialize(Npix, R);
params.cfg = cfg;
params.vectors = vectors;
params.split = 1;
it_max = 500;
xld_input = cat(3,xld_H,xld_V);
mask = mx.^2 + my.^2 + mz.^2 > 0;

clear xld_H xld_V
clear mx my mz
clear circ height
close all

%% 
for ii=1:20
    mx_init = (rand(Npix)-0.5)*2;
    my_init = (rand(Npix)-0.5)*2;
    mz_init = (rand(Npix)-0.5)*2;
    
    tic
    [mx_out,my_out,mz_out,E] = xld_recons_all_init(xld_input,R,10,mask,params,pol_vec,mx_init,my_init,mz_init);
    toc

    filename = sprintf('%s/rand%02d.mat',path_name,ii);

    save(filename,'mx_out','my_out','mz_out','E','cfg','vectors')  
end
