function [error_image, theta_error, mag_error_image] = config_error(mx_rec,my_rec,mz_rec,mx_init,my_init,mz_init)
    recons_1 = cat(4, mx_rec, my_rec, mz_rec);
    recons_2 = -recons_1;
    
    config = cat(4, mx_init, my_init, mz_init);
    error_1 = sum((recons_1 - config).^2,4);
    error_2 = sum((recons_2 - config).^2,4);
    
    error_image = error_1;
    error_image(error_1 > error_2) = error_2(error_1 > error_2);
    
    mag_c = sqrt(sum(config.^2,4));
    mag_rc = sqrt(sum(recons_1.^2,4));
    mag_error_image = abs(mag_c - mag_rc);
    
    dot_prod = abs(dot(recons_1,config,4));
    div = mag_c.* mag_rc;
    argument = dot_prod ./ div;
    theta_error = real(acosd(argument));
    
    % error in each component
    % angle to tilt axis?? 
end
