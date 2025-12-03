%% 1.0 Load Data
load('off_resonance_projections.mat');
load('aligned_proj_tilt1.mat');
load('aligned_proj_tilt2.mat');
load('aligned_proj_tilt3.mat');
load('aligned_proj_tilt4.mat');

%% 2.0 Reconstruct electron density tomogram and use as mask
Npix = [650, 650, 760]; % size of the 3D reconstruction volume

[cfg, vectors] = astra_initialize(Npix,aligned_q.R_ch);
phase_reconstruction = FBP(aligned_q.data_ch_shift, cfg, vectors);
charge_reeconstuction = -phase_reconstruction*factor_edensity;

mask = and(phase_reconstruction < -0.015, phase_reconstruction > -0.01825);
clear tomo_q_full

%% 3.0 Reconstruction Parameters
pol_H = 0;
pol_V = 90;
V_purity = 0.64;
n_recons = 5;

R_r = aligned_r.R_lv;
data_r = aligned_r.data_lh_shift-aligned_r.data_lv_shift;

R_a = aligned_a.R_lv;
data_a = aligned_a.data_lh_shift-aligned_a.data_lv_shift;

R_b = aligned_b.R_lv;
data_b = aligned_b.data_lh_shift-aligned_b.data_lv_shift;

R_c = aligned_c.R_lv;
data_c = aligned_c.data_lh_shift-aligned_c.data_lv_shift;

R_all = cat(3,R_r,R_a,R_b,R_c);
data_all = cat(3,data_r,data_a,data_b,data_c);


% removed unused variables from RAM
clear aligned_r aligned_a aligned_b aligned_c R_r R_a R_b R_c
clear data_r data_a data_b data_c

%% 3.1 Full XLD reconstruction
[cfg, vectors] = astra_initialize(Npix,R_all);
params.cfg = cfg;
params.vectors = vectors;
params.split = 2;

for ii=1:n_recons
    init_x = single(0.1.*rand(Npix)-0.05);
    init_y = single(0.1.*rand(Npix)-0.05);
    init_z = single(0.1.*rand(Npix)-0.05);
    it_max = 50;

    tic
    [mx,my,mz,E] = xld_recons_diff_init_state(data_all,R_all,it_max,mask,params,pol_H,pol_V,V_purity,init_x,init_y,init_z);
    toc
    mx = single(0.1.*ones(Npix)-0.05);
    my = single(0.1.*ones(Npix)-0.05);
    mz = single(0.1.*ones(Npix)-0.05);
    
    filename = sprintf('v2o5_xld%02d.mat',ii);
    save(filename,'mx','my','mz','E','cfg','vectors','-v7.3')
    clear mx my mz
end

%% 3.2 Load and average reconstruction

recons_stack = zeros(n_recons,3,330,330,540);
for ii=1:n_recons
    filename = sprintf('v2o5_xld%02d.mat',ii);
    load(filename,'mx','my','mz')
    recons_stack(ii,1,:,:,:) = mx(151:480,151:480,171:710);
    recons_stack(ii,2,:,:,:) = my(151:480,151:480,171:710);
    recons_stack(ii,3,:,:,:) = mz(151:480,151:480,171:710);
end

[avg_x, avg_y, avg_z] = average_recons(recons_stack);
save('v2o5_average.mat','mx','my','mz','-v7.3')