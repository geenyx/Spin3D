function [final_x, final_y, final_z] = average_recons(recons_stack)
n_recons = size(recons_stack,1);
for ii=1:n_recons
    %%%%%%%%%%%%%% REMOVE ZERO PADDING
    temp_x = recons_stack(ii,1,:,:,:);
    temp_y = recons_stack(ii,2,:,:,:);
    temp_z = recons_stack(ii,3,:,:,:);
    
    %%%%%%%%%%%%%% X ALIGNMENT
    indices = find(temp_x<0);
    temp_x(indices) = -temp_x(indices);
    temp_y(indices) = -temp_y(indices);
    temp_z(indices) = -temp_z(indices);
       
    recons_stack_xalign(ii,1,:,:,:) = temp_x;
    recons_stack_xalign(ii,2,:,:,:) = temp_y;
    recons_stack_xalign(ii,3,:,:,:) = temp_z;
    
    %%%%%%%%%%%%%% Y ALIGNMENT
    indices = find(temp_y<0);
    temp_x(indices) = -temp_x(indices);
    temp_y(indices) = -temp_y(indices);
    temp_z(indices) = -temp_z(indices);
    
    recons_stack_yalign(ii,1,:,:,:) = temp_x;
    recons_stack_yalign(ii,2,:,:,:) = temp_y;
    recons_stack_yalign(ii,3,:,:,:) = temp_z;
    
    %%%%%%%%%%%%%% Z ALIGNMENT
    indices = find(temp_z<0);
    temp_x(indices) = -temp_x(indices);
    temp_y(indices) = -temp_y(indices);
    temp_z(indices) = -temp_z(indices);
    
    recons_stack_zalign(ii,1,:,:,:) = temp_x;
    recons_stack_zalign(ii,2,:,:,:) = temp_y;
    recons_stack_zalign(ii,3,:,:,:) = temp_z;
end

%%%%%%%%% X ALIGNED MEAN 
mean_x_xalign = squeeze(mean(recons_stack_xalign(:,1,:,:,:),1));
mean_y_xalign = squeeze(mean(recons_stack_xalign(:,2,:,:,:),1));
mean_z_xalign = squeeze(mean(recons_stack_xalign(:,3,:,:,:),1));

std_x_xalign = squeeze(std(recons_stack_xalign(:,1,:,:,:),1));
std_y_xalign = squeeze(std(recons_stack_xalign(:,2,:,:,:),1));
std_z_xalign = squeeze(std(recons_stack_xalign(:,3,:,:,:),1));

%%%%%%%%% Y ALIGNED MEAN 
mean_x_yalign = squeeze(mean(recons_stack_yalign(:,1,:,:,:),1));
mean_y_yalign = squeeze(mean(recons_stack_yalign(:,2,:,:,:),1));
mean_z_yalign = squeeze(mean(recons_stack_yalign(:,3,:,:,:),1));

std_x_yalign = squeeze(std(recons_stack_yalign(:,1,:,:,:),1));
std_y_yalign = squeeze(std(recons_stack_yalign(:,2,:,:,:),1));
std_z_yalign = squeeze(std(recons_stack_yalign(:,3,:,:,:),1));

%%%%%%%%% Z ALIGNED MEAN 
mean_x_zalign = squeeze(mean(recons_stack_zalign(:,1,:,:,:),1));
mean_y_zalign = squeeze(mean(recons_stack_zalign(:,2,:,:,:),1));
mean_z_zalign = squeeze(mean(recons_stack_zalign(:,3,:,:,:),1));

std_x_zalign = squeeze(std(recons_stack_zalign(:,1,:,:,:),1));
std_y_zalign = squeeze(std(recons_stack_zalign(:,2,:,:,:),1));
std_z_zalign = squeeze(std(recons_stack_zalign(:,3,:,:,:),1));

%%%%%%%%% X ALIGNED ERRORS 
theta_error_stack_xaligned = zeros(n_recons,330,330,540);
for ii=1:n_recons
    mx_temp = squeeze(recons_stack_xalign(ii,1,:,:,:));
    my_temp = squeeze(recons_stack_xalign(ii,2,:,:,:));
    mz_temp = squeeze(recons_stack_xalign(ii,3,:,:,:));

    [~, theta_error, ~] = config_error(mx_temp,my_temp,mz_temp,mean_x_xalign,mean_y_xalign,mean_z_xalign);
    theta_error_stack_xaligned(ii,:,:,:) = theta_error;
end

stdev_theta_error_xalign = sqrt(squeeze(sum(theta_error_stack_xaligned.^2,1))./n_recons);
clear theta_error_stack_xaligned mx_temp my_temp mz_temp

%%%%%%%%% Y ALIGNED ERRORS 
theta_error_stack_yaligned = zeros(n_recons,330,330,540);
for ii=1:n_recons
    mx_temp = squeeze(recons_stack_yalign(ii,1,:,:,:));
    my_temp = squeeze(recons_stack_yalign(ii,2,:,:,:));
    mz_temp = squeeze(recons_stack_yalign(ii,3,:,:,:));

    [~, theta_error, ~] = config_error(mx_temp,my_temp,mz_temp,mean_x_yalign,mean_y_yalign,mean_z_yalign);
    theta_error_stack_yaligned(ii,:,:,:) = theta_error;
end
stdev_theta_error_yalign = sqrt(squeeze(sum(theta_error_stack_yaligned.^2,1))./n_recons);
clear theta_error_stack_yaligned mx_temp my_temp mz_temp

%%%%%%%%% Z ALIGNED ERRORS 
theta_error_stack_zaligned = zeros(n_recons,330,330,540);
for ii=1:n_recons
    mx_temp = squeeze(recons_stack_zalign(ii,1,:,:,:));
    my_temp = squeeze(recons_stack_zalign(ii,2,:,:,:));
    mz_temp = squeeze(recons_stack_zalign(ii,3,:,:,:));

    [~, theta_error, ~] = config_error(mx_temp,my_temp,mz_temp,mean_x_zalign,mean_y_zalign,mean_z_zalign);
    theta_error_stack_zaligned(ii,:,:,:) = theta_error;
end
stdev_theta_error_zalign = sqrt(squeeze(sum(theta_error_stack_zaligned.^2,1))./n_recons);
clear theta_error_stack_zaligned mx_temp my_temp mz_temp

final_x = mean_x_xalign;
final_y = mean_y_xalign;
final_z = mean_z_xalign;

% start with x, and then replace the values where yalign has smaller
% deviation, then repeat for z
yltx = stdev_theta_error_yalign < stdev_theta_error_xalign;
zlty = and(stdev_theta_error_zalign < stdev_theta_error_yalign,stdev_theta_error_zalign < stdev_theta_error_xalign);

% replace x with y
final_x(yltx) = mean_x_yalign(yltx);
final_y(yltx) = mean_y_yalign(yltx);
final_z(yltx) = mean_z_yalign(yltx);

% replace y with z
final_x(zlty) = mean_x_zalign(zlty);
final_y(zlty) = mean_y_zalign(zlty);
final_z(zlty) = mean_z_zalign(zlty);
