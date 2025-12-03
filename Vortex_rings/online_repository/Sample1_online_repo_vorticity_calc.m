%% Preparing files for paraview
close all
clear 

load('Sample1_magnetic_recons_mag_3comp.mat');

load('Sample1_nonmagnetic_recons.mat');

%% Crop to remove zeros
mx_out(:,[1:70 418:487],:)=[];
my_out(:,[1:70 418:487],:)=[];
mz_out(:,[1:70 418:487],:)=[];
mx_out([1:70 418:487],:,:)=[];
my_out([1:70 418:487],:,:)=[];
mz_out([1:70 418:487],:,:)=[];

%% for paraview:
freq_scale = 0.07;
freq_scale_z = 0.07;
mx_ = filter_x(mx_out,freq_scale);
mx_ = filter_y(mx_,freq_scale);
mx_ = filter_z(mx_,freq_scale);
my_ = filter_x(my_out,freq_scale);
my_ = filter_y(my_,freq_scale);
my_ = filter_z(my_,freq_scale);
mz_ = filter_x(mz_out,freq_scale_z);
mz_ = filter_y(mz_,freq_scale_z);
mz_ = filter_z(mz_,freq_scale_z);



%%

slicenum = [1:440];
mx = mx_(:,:,slicenum);
my = my_(:,:,slicenum);
mz = mz_(:,:,slicenum);
tomogram_beta_ = tomogram_beta(:,:,slicenum);

mx_small = imresize(mx,1);
my_small = imresize(my,1);
mz_small = imresize(mz,1);
tomo_recons = imresize(tomogram_beta_,1);

X = [1:size(mx_small,1)];
Y = [1:size(mx_small,2)];
Z = [1:size(mx_small,3)];
[x,y,z] = meshgrid(X,Y,Z);



%% Mask using the non-magnetic tomogram:
 mask_mat = abs(tomo_recons)>3e-6;

 tomo_recons = tomo_recons.*mask_mat;
 mx_small = mx_small.*mask_mat;
 my_small = my_small.*mask_mat;
 mz_small = mz_small.*mask_mat;
 divM = divM.*mask_mat;
 
 
 %% Calculate magnetic vorticity:
 
 m_abs = sqrt(mx_small.^2+my_small.^2+mz_small.^2);
 
 % Define your vector:
 
 nx = mx_small./m_abs;
 ny = my_small./m_abs;
 nz = mz_small./m_abs;
 
 % First calculate the different gradients:
 
 
 [nx_x, nx_y, nx_z] = gradient(nx);
 [ny_x, ny_y, ny_z] = gradient(ny);
 [nz_x, nz_y, nz_z] = gradient(nz);
 
 mag_vort_x = 2*nx.*(ny_y.*nz_z - ny_z.*nz_y) + ...
            2*ny.*(nz_y.*nx_z - nz_z.*nx_y) + ...
            2*nz.*(nx_y.*ny_z - nx_z.*ny_y);
 
 mag_vort_y = 2*nx.*(ny_z.*nz_x - ny_x.*nz_z) + ...
            2*ny.*(nz_z.*nx_x - nz_x.*nx_z) + ...
            2*nz.*(nx_z.*ny_x - nx_x.*ny_z);
        
 mag_vort_z = 2*nx.*(ny_x.*nz_y - ny_y.*nz_x) + ...
            2*ny.*(nz_x.*nx_y - nz_y.*nx_x) + ...
            2*nz.*(nx_x.*ny_y - nx_y.*ny_x);
 
 n_abs = sqrt(nx.^2+ny.^2+nz.^2);
 %% Have a look at the vorticity:
 
 mag_vort_x(isnan(mag_vort_x))=0;
 mag_vort_y(isnan(mag_vort_y))=0;
 mag_vort_z(isnan(mag_vort_z))=0;
 n_abs(isnan(n_abs))=0;
 vort_abs = sqrt(mag_vort_x.^2+mag_vort_y.^2+mag_vort_z.^2);
 
 vort_axis = [-0.05 0.05];
 
 for ii = 100
 figure(1)
 subplot(3,1,1)
 imagesc(mag_vort_x(:,:,ii))
 axis xy equal tight
 colorbar
 caxis(vort_axis)
 subplot(3,1,2)
 imagesc(mag_vort_y(:,:,ii))
 axis xy equal tight
 colorbar
 caxis(vort_axis)
 subplot(3,1,3)
 imagesc(mag_vort_z(:,:,ii))
 axis xy equal tight
 colorbar
 caxis(vort_axis)

 end
 

 
 %% Save vtk files for paraview
 
 vtkwrite('Sample_1_remanent_asgrown_vorticity.vtk','structured_grid',x,y,z,'vectors','mag_vorticity',mag_vort_x,mag_vort_y,mag_vort_z,'scalars','m_abs',n_abs);
 vtkwrite('Sample_1_remanent_asgrown_magnetisation_and_vorticty.vtk','structured_grid',x,y,z,'vectors','magnetisation',mx_small,my_small,mz_small,'vectors','mag_vorticity',mag_vort_x,mag_vort_y,mag_vort_z,'scalars','m_abs',n_abs);
 vtkwrite('Sample_1_remanent_asgrown_magnetisation.vtk','structured_grid',x,y,z,'vectors','magnetisation',mx_small,my_small,mz_small,'scalars','m_abs',n_abs);
