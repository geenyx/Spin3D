%   tomo_filtered = filter_y(tomo,freq_scale)
%   Receives a tomogram and applies filtering only along the third index of
%   a 3D matrix. 
% Inputs:
%   tomo        Input tomogram matrix
%   freq_scale  Frequency cutoff
% Manuel Guizar Feb 14, 2016

function tomo_filtered = filter_x(tomo_filtered,freq_scale)


%%% Filter a tomogram along the last index %%%
% freq_scale = 0.1;
% 
% %%% Create example %%%
% N = 200;
% tomo_filtered = zeros([N N N]);
% for ii = 1:N
%     tomo_filtered(ii,:,:) = phantom(N);
% end
% tomo_filtered = repmat(tomo_filtered,[1 1 N]);
%%%%%%%%%%%%%%%%%%%%%%

%%% Create filter %%%

Nfilt = size(tomo_filtered,2);

filt = zeros([Nfilt 1]);
d = freq_scale;

w = [-Nfilt/2:Nfilt/2-1]+mod(Nfilt/2,2);
w = 2*pi*w/Nfilt;

filt = (1+cos(w/d)) / 2;
filt(abs(w)/d>pi) = 0;
filt = ifftshift(filt);

% figure(1);
% plot(filt);
%%%%%%
%

filt3D = repmat(reshape(filt,[1,Nfilt,1]),[size(tomo_filtered,1) 1 size(tomo_filtered,3)]);

tomo_filtered = real(ifft( fft(tomo_filtered,[],2).*filt3D ,[],2));



% %% Image %%%
% figure(1)
% imagesc(squeeze(tomo_filtered(2,:,:)));
% % imagesc(squeeze(filt3D(1,:,:)));
% colorbar

