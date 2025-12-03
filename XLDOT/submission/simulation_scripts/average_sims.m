%% XLD reconstruction
clear
close all

addpath functions

%% Load dataset, savepath
path_name = 'output';
regimes = {'smoothly_varying','granular'};
methods = {'single_axis','dual_axis','triple_axis'};

for regime_idx=1:length(regimes)
    clearvars -except path_name regimes methods regime_idx 
for method_idx=1:length(methods) 
%% Reference
load(sprintf('reference_structures/%s.mat',regimes{regime_idx}));

%% Averaging
for n_recons=[20] % change to [5,10,15,20] to calculate different number of averages
recons_stack_xalign = zeros(n_recons,3,size(mx,1),size(mx,2),size(mx,3));
recons_stack_yalign = zeros(n_recons,3,size(mx,1),size(mx,2),size(mx,3));
recons_stack_zalign = zeros(n_recons,3,size(mx,1),size(mx,2),size(mx,3));

for ii=1:n_recons
    filename = sprintf('%s/%s/%s/rand%02d.mat',path_name,regimes{regime_idx},methods{method_idx},ii);
    load(filename,'mx_out','my_out','mz_out','E','cfg','vectors')
    stack(ii,1,:,:,:) = mx_out;
    stack(ii,2,:,:,:) = my_out;
    stack(ii,3,:,:,:) = mz_out;

    temp_x = mx_out;
    temp_y = my_out;
    temp_z = mz_out;
    mask = temp_x.^2 + temp_y.^2 + temp_z.^2 > 0;
    %%%%%%%%%%%%%% X ALIGNMENT
    indices = find(temp_x<0);
    temp_x(indices) = -temp_x(indices);
    temp_y(indices) = -temp_y(indices);
    temp_z(indices) = -temp_z(indices);
    
    temp_x_filt = temp_x;
    temp_y_filt = temp_y;
    temp_z_filt = temp_z;
            
    boundary_mask = and(temp_x_filt == 0, mask ==1);
    temp_x_filt(boundary_mask) = temp_x(boundary_mask);
    boundary_mask = and(temp_y_filt == 0, mask ==1);
    temp_y_filt(boundary_mask) = temp_x(boundary_mask);
    boundary_mask = and(temp_z_filt == 0, mask ==1);
    temp_z_filt(boundary_mask) = temp_x(boundary_mask);
    
    recons_stack_xalign(ii,1,:,:,:) = temp_x_filt;
    recons_stack_xalign(ii,2,:,:,:) = temp_y_filt;
    recons_stack_xalign(ii,3,:,:,:) = temp_z_filt;
    
    %%%%%%%%%%%%%% Y ALIGNMENT
    indices = find(temp_y<0);
    temp_x(indices) = -temp_x(indices);
    temp_y(indices) = -temp_y(indices);
    temp_z(indices) = -temp_z(indices);
    

    temp_x_filt = temp_x;
    temp_y_filt = temp_y;
    temp_z_filt = temp_z;
    
    boundary_mask = and(temp_x_filt == 0, mask ==1);
    temp_x_filt(boundary_mask) = temp_x(boundary_mask);
    boundary_mask = and(temp_y_filt == 0, mask ==1);
    temp_y_filt(boundary_mask) = temp_x(boundary_mask);
    boundary_mask = and(temp_z_filt == 0, mask ==1);
    temp_z_filt(boundary_mask) = temp_x(boundary_mask);
    
    recons_stack_yalign(ii,1,:,:,:) = temp_x_filt;
    recons_stack_yalign(ii,2,:,:,:) = temp_y_filt;
    recons_stack_yalign(ii,3,:,:,:) = temp_z_filt;
    
    %%%%%%%%%%%%%% Z ALIGNMENT
    indices = find(temp_z<0);
    temp_x(indices) = -temp_x(indices);
    temp_y(indices) = -temp_y(indices);
    temp_z(indices) = -temp_z(indices);
    
    temp_x_filt = temp_x;
    temp_y_filt = temp_y;
    temp_z_filt = temp_z;
    
    boundary_mask = and(temp_x_filt == 0, mask ==1);
    temp_x_filt(boundary_mask) = temp_x(boundary_mask);
    boundary_mask = and(temp_y_filt == 0, mask ==1);
    temp_y_filt(boundary_mask) = temp_x(boundary_mask);
    boundary_mask = and(temp_z_filt == 0, mask ==1);
    temp_z_filt(boundary_mask) = temp_x(boundary_mask);
    
    recons_stack_zalign(ii,1,:,:,:) = temp_x_filt;
    recons_stack_zalign(ii,2,:,:,:) = temp_y_filt;
    recons_stack_zalign(ii,3,:,:,:) = temp_z_filt;
end

clear temp_x temp_y temp_z indices temp_x_filt temp_y_filt temp_z_filt boundary_mask

%%%%%%%%% X ALIGNED median 
median_x_xalign = squeeze(median(recons_stack_xalign(:,1,:,:,:),1));
median_y_xalign = squeeze(median(recons_stack_xalign(:,2,:,:,:),1));
median_z_xalign = squeeze(median(recons_stack_xalign(:,3,:,:,:),1));

std_x_xalign = squeeze(std(recons_stack_xalign(:,1,:,:,:),1));
std_y_xalign = squeeze(std(recons_stack_xalign(:,2,:,:,:),1));
std_z_xalign = squeeze(std(recons_stack_xalign(:,3,:,:,:),1));

%%%%%%%%% Y ALIGNED median 
median_x_yalign = squeeze(median(recons_stack_yalign(:,1,:,:,:),1));
median_y_yalign = squeeze(median(recons_stack_yalign(:,2,:,:,:),1));
median_z_yalign = squeeze(median(recons_stack_yalign(:,3,:,:,:),1));

std_x_yalign = squeeze(std(recons_stack_yalign(:,1,:,:,:),1));
std_y_yalign = squeeze(std(recons_stack_yalign(:,2,:,:,:),1));
std_z_yalign = squeeze(std(recons_stack_yalign(:,3,:,:,:),1));

%%%%%%%%% Z ALIGNED median 
median_x_zalign = squeeze(median(recons_stack_zalign(:,1,:,:,:),1));
median_y_zalign = squeeze(median(recons_stack_zalign(:,2,:,:,:),1));
median_z_zalign = squeeze(median(recons_stack_zalign(:,3,:,:,:),1));

std_x_zalign = squeeze(std(recons_stack_zalign(:,1,:,:,:),1));
std_y_zalign = squeeze(std(recons_stack_zalign(:,2,:,:,:),1));
std_z_zalign = squeeze(std(recons_stack_zalign(:,3,:,:,:),1));

%%%%%%%%% X ALIGNED ERRORS 
theta_error_stack_xaligned = zeros(n_recons,size(mx_out,1),size(mx_out,2),size(mx_out,3));
for ii=1:n_recons
    fprintf(sprintf('Calculating errors (x aligned): %d\n',ii))
    mx_temp = squeeze(recons_stack_xalign(ii,1,:,:,:));
    my_temp = squeeze(recons_stack_xalign(ii,2,:,:,:));
    mz_temp = squeeze(recons_stack_xalign(ii,3,:,:,:));

    [~, theta_error, ~] = config_error(mx_temp,my_temp,mz_temp,median_x_xalign,median_y_xalign,median_z_xalign);
    theta_error_stack_xaligned(ii,:,:,:) = theta_error;
end
median_theta_error_xalign = squeeze(median(theta_error_stack_xaligned,1));
stdev_theta_error_xalign = sqrt(squeeze(sum(theta_error_stack_xaligned.^2,1))./n_recons);
clear theta_error_stack_xaligned mx_temp my_temp mz_temp

%%%%%%%%% Y ALIGNED ERRORS 
theta_error_stack_yaligned = zeros(n_recons,size(mx_out,1),size(mx_out,2),size(mx_out,3));
for ii=1:n_recons
    fprintf(sprintf('Calculating errors (y aligned): %d\n',ii))
    mx_temp = squeeze(recons_stack_yalign(ii,1,:,:,:));
    my_temp = squeeze(recons_stack_yalign(ii,2,:,:,:));
    mz_temp = squeeze(recons_stack_yalign(ii,3,:,:,:));

    [~, theta_error, ~] = config_error(mx_temp,my_temp,mz_temp,median_x_yalign,median_y_yalign,median_z_yalign);
    theta_error_stack_yaligned(ii,:,:,:) = theta_error;
end
median_theta_error_yalign = squeeze(median(theta_error_stack_yaligned,1));
stdev_theta_error_yalign = sqrt(squeeze(sum(theta_error_stack_yaligned.^2,1))./n_recons);
clear theta_error_stack_yaligned mx_temp my_temp mz_temp

%%%%%%%%% Z ALIGNED ERRORS 
theta_error_stack_zaligned = zeros(n_recons,size(mx_out,1),size(mx_out,2),size(mx_out,3));
for ii=1:n_recons
    fprintf(sprintf('Calculating errors (z aligned): %d\n',ii))
    mx_temp = squeeze(recons_stack_zalign(ii,1,:,:,:));
    my_temp = squeeze(recons_stack_zalign(ii,2,:,:,:));
    mz_temp = squeeze(recons_stack_zalign(ii,3,:,:,:));

    [~, theta_error, ~] = config_error(mx_temp,my_temp,mz_temp,median_x_zalign,median_y_zalign,median_z_zalign);
    theta_error_stack_zaligned(ii,:,:,:) = theta_error;
end
median_theta_error_zalign = squeeze(median(theta_error_stack_zaligned,1));
stdev_theta_error_zalign = sqrt(squeeze(sum(theta_error_stack_zaligned.^2,1))./n_recons);
clear theta_error_stack_zaligned mx_temp my_temp mz_temp

final_x = median_x_xalign;
final_y = median_y_xalign;
final_z = median_z_xalign;
final_std = stdev_theta_error_xalign;
% start with x, and then replace the values where yalign has smaller
% deviation, then repeat for z

%%%%%%%%%%%
%%%%%%%%%%%
yltx = stdev_theta_error_yalign < stdev_theta_error_xalign;
zlty = and(stdev_theta_error_zalign < stdev_theta_error_yalign,stdev_theta_error_zalign < stdev_theta_error_xalign);
% replace x with y
final_x(yltx) = median_x_yalign(yltx);
final_y(yltx) = median_y_yalign(yltx);
final_z(yltx) = median_z_yalign(yltx);
final_std(yltx) = stdev_theta_error_yalign(yltx);
% replace y with z
final_x(zlty) = median_x_zalign(zlty);
final_y(zlty) = median_y_zalign(zlty);
final_z(zlty) = median_z_zalign(zlty);
final_std(zlty) = stdev_theta_error_zalign(zlty);

%% statistics from averaging
error_stack =  zeros(n_recons,size(mx,1),size(mx,2),size(mx,3));

for jj=1:n_recons
    rec_mx = squeeze(stack(jj,1,:,:,:));
    rec_my = squeeze(stack(jj,2,:,:,:));
    rec_mz = squeeze(stack(jj,3,:,:,:));

    [~,rec_theta_error,~] = config_error(rec_mx,rec_my,rec_mz,final_x,final_y,final_z);
    error_stack(jj,:,:,:) = rec_theta_error;
end

uncertainty = squeeze(mean(error_stack,1));

%% configuration error
[error,theta_error, ~]= config_error(final_x,final_y,final_z,mx,my,mz);

%% save
s_filename = sprintf('%s/%s/%s/average%02d.mat',path_name,regimes{regime_idx},methods{method_idx},n_recons);

save(s_filename, 'final_x', 'final_y', 'final_z','error','theta_error','uncertainty') 
end
end
end
