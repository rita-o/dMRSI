clear all
clc
close all

toolbox_path='/home/localadmin/Documents/Rita/Toolboxes';
addpath(genpath(fullfile(toolbox_path,'Tensor-MP-PCA-main')))
addpath(genpath(fullfile(toolbox_path,'MP-PCA-main')))

path_save = '/home/localadmin/Documents/Rita/Data/dMRSI_denoising';

%% Choose dataset

files = dir(fullfile(path_save,'*.mat'));

%for f=1:length(files)
%f=11
%load(fullfile(path_save,files(f).name))
load(fullfile(path_save,'DSI94R110_1139.mat'))

data_3D = fftshift(fftshift(kspaceData,2),3);
% figure; 
% imagesc(squeeze(real(k(10,:,:)))); 
% colormap('gray')
% 
% imgFIDs = ifft(ifft(k,[],2),[],3);
% data_3D = ifftshift(ifftshift(imgFIDs,2),3);    
% 
% figure; 
% imagesc(squeeze(real(data_3D(10,:,:)))); 
% colormap('gray')


%end


% get info
npoints  = 768;         % number of FID points (taken from methods file)
B0 = 9.4;               % Tesla
gamma = 42.57747892e6;  % Hz/T for 1H
f0 =B0*gamma*1e-6;      % central frequency in MHz
swh = 5000;              % Hz (taken from methods file)
sw = swh*2;
dw = 1/sw;              % dwell time (in s)
grpdly   = 76;          % group delay (taken from acq file)

fmax      = swh/2;
%fmax = swh;
freq_axis = [fmax:-2*fmax/(npoints-1):-fmax];
ppm_axis  = freq_axis/f0+4.7;

%%
% example dataset:
fid = squeeze(data_3D(:,13,13)).';

% appodization
% time_axis = 0:dw:(npoints-1)*dw;
% offset_hz = -2 * f0 ;
% fid = fid .* exp(1i* offset_hz * 2 * pi *time_axis);                   
%fid = fid - mean(fid(end-64:end));         
% 


%group delay
if grpdly>0 
    fid = [fid(grpdly+1:end), zeros(1,grpdly)];
end
% 
fid = conj(fid); 

spec = fftshift(fft(fid, npoints, 2), 2);


% win = (ppm_axis>1) & (ppm_axis<4);
% phi0 = -angle(sum(spec(win)));
% spec = spec * exp(1i*phi0);

% figure
% %subplot(1,2,1)
% plot(real(fid)) 
% 
%subplot(1,2,2)
figure
plot(ppm_axis,real(spec)) 
xlabel('ppm')
set(gca, 'XDir','reverse')



%% Denoise

data_3D_perm = permute(data_3D, [2 3 1]);
[data_3D_dn_tmp, sigma, P] = denoise(data_3D_perm,[5,5]);
data_3D_dn = permute(data_3D_dn_tmp, [3 1 2]);

data_3D_perm = permute(data_3D, [2 3 1]);
[data_3D_dn_tmp, sigma, P] = denoise_recursive_tensor(data_3D_perm,[5,5]);
data_3D_dn_t = permute(data_3D_dn_tmp, [3 1 2]);

%%
[spec]                 = prepare_spectrum(data_3D, grpdly, dw, f0);
[spec_dn]              = prepare_spectrum(data_3D_dn, grpdly, dw, f0);
[spec_dn_t]              = prepare_spectrum(data_3D_dn_t, grpdly, dw, f0);


make_figure(spec,spec_dn,spec_dn_t, ppm_axis)
saveas(gcf,fullfile(path_save,'exp.png'))

% figure
% plot(ppm_axis,real(spec(:,1,1))/1e6,'k')
% hold on% plot(ppm_axis,real(spec_dn(:,1,1))/1e6,'b')
% xlabel('ppm')
% set(gca, 'XDir','reverse')

% make_figure(ft_sum, ft_dn_tmp2D_sum,ft_dn_tmp4D_sum, ppmscale, bvals, deltas,...
%     {'Original','t-MP-PCA 2D','res t-MP-PCA 2D','t-MP-PCA 4D','res t-MP-PCA 4D'},...
%     fullfile(path_save,'Original_vs_tMPPCA2D_vs_tMPPCA4D_summed'))


%% Auxiliar Function

% Prepare spectrum
function [spec] = prepare_spectrum(data, grpdly, dw, f0)

npoints  = size(data,1);
nx       = size(data,2);
ny       = size(data,3);
for x = 1:nx
    for y = 1:ny
        fid = data(:,x,y).';

        % appodization
        % time_axis = 0:dw:(npoints-1)*dw;
        % offset_hz = -2 * f0 ;
        % fid = fid .* exp(1i* offset_hz * 2 * pi *time_axis);                   
        fid = conj(fid); 
        % fid = fid - mean(fid(end-64:end));            % or: fid = fid - mean(fid);
      
        % group delay
        if grpdly>0 
            fid = [fid(grpdly+1:end), zeros(1,grpdly)];
        end
        
        %spec = fftshift(fft(fid, [], 2), 2);
        % win = (ppm_axis>1) & (ppm_axis<4);
        % phi0 = -angle(sum(spec(win)));
        % spec = spec * exp(1i*phi0);

        spec(:,x,y) =fftshift(fft(fid,npoints,2),2);
    end
end

end

% Make figure sum
function [] = make_figure(data1, data2, data3, ppmscale)

close all

% Get residuals
res2 = real(data2-data1);
res3 = real(data3-data1);

% Define dimensionss
nx       = size(data1,2);
ny       = size(data2,3);
nb_plots = 5;
rows = round(linspace(1,nx,nb_plots));
cols = round(linspace(1,ny,nb_plots));
[RR, CC] = ndgrid(rows, cols);
list_pairs = [CC(:) RR(:)];  

figure;
tiledlayout(nb_plots, nb_plots, 'Padding', 'compact', 'TileSpacing', 'compact');
set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 40, 20]);

for i=1:nb_plots*nb_plots
    pair = list_pairs(i,:);
    nexttile
    plot(ppmscale,real(data1(:,pair(1),pair(2)))/1e4,'k')
    hold on
    plot(ppmscale,real(data2(:,pair(1),pair(2)))/1e4,'b')
    plot(ppmscale,real(res2(:,pair(1),pair(2)))/1e4-2.5,'b')

    plot(ppmscale,real(data3(:,pair(1),pair(2)))/1e4,'r')
    plot(ppmscale,real(res3(:,pair(1),pair(2)))/1e4-4.5,'r')

    xlabel('ppm')
    ylabel('Signal×10^4')
    ylim([-6 4])
    %ylim([-200 100])
    xlim([0,6])
    title(sprintf('Image position %d,%d',pair(1),pair(2)))
    if i==nb_plots*nb_plots
        legend({'Original','MP-PCA','MP-PCA res','t-MP-PCA','t-MP-PCA res'},'Location','northwest')
    end
    set(gca, 'XDir','reverse')
    grid on
end

end
