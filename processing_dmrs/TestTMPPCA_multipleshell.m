clear all
clc
close all

toolbox_path='/home/localadmin/Documents/Rita/Toolboxes';
addpath(genpath(fullfile(toolbox_path,'Tensor-MP-PCA-main')))
addpath(genpath(fullfile(toolbox_path,'MP-PCA-main')))

path_save = '/home/localadmin/Documents/Rita/Data/dMRI_dMRS_Pilot_20250424/derivatives/preprocessed/sub-01/ses-01/dmrs';

%% Choose dataset

path_data = '/home/localadmin/Documents/Rita/Data/dMRI_dMRS_Pilot_20250424/derivatives/preprocessed/sub-01/ses-01/dmrs/combined';
load(fullfile(path_data,'sub-01_ses-01_combined_dmrs_allDelta-allb_FID.mat'))
load(fullfile(path_data,'sub-01_ses-01_combined_dmrs_allDelta-allb_info.mat'))

% sizes
nshells = size(data,1);
ndeltas = size(data,2);
nshots  = size(data,3);
nt      = repmat(nshots, 1, nshells)';  % 1 × nshells
npoints = size(data,4);

% get info
B0field = data_info.centralFreq*1e-6;
sw      = 1/data_info.dwelltime;
fmax     = sw/2;
f        = [fmax:-2*fmax/(npoints-1):-fmax];
ppmscale = f/B0field+4.7;
ppm_mask = ppmscale <= 10.66 & ppmscale >= 8.2;
grpdly   = 76 % check this
grpdly   = repmat(grpdly, 1, nshells);  % 1 × nshells
bvals    = data_info.bvals;
deltas   = data_info.deltas;

% get data (4D)
data_4D = data;

% concatenate data (2D)
[nshells, ndeltas, nshots, npoints] = size(data_4D);
data_2D = [];
for d = 1:ndeltas
    for b=1:nshells
        for s=1:nshots
                data_2D = [data_2D; squeeze(data_4D(b, d, s, :)).'];
        end
    end
end

% get data of largest and shortest diff time data 
data_4D_DeltaMax = data_4D(:,end,:,:);
data_4D_DeltaMin = data_4D(:,1,:,:);

% get data of largest and shortest b value data 
data_4D_b1 = data_4D(1,:,:,:);
data_4D_b2 = data_4D(2,:,:,:);
data_4D_b3 = data_4D(3,:,:,:);
data_4D_b4 = data_4D(4,:,:,:);
data_4D_b5 = data_4D(end,:,:,:);

% concatenate data (2D)
[nshells, ndeltas, nshots, npoints] = size(data_4D_b5);

data_2D_b1 = [];
for d = 1:ndeltas
    for b=1:nshells
        for s=1:nshots
                data_2D_b1 = [data_2D_b1; squeeze(data_4D_b1(b, d, s, :)).'];
        end
    end
end
data_2D_b2 = [];
for d = 1:ndeltas
    for b=1:nshells
        for s=1:nshots
                data_2D_b2 = [data_2D_b2; squeeze(data_4D_b2(b, d, s, :)).'];
        end
    end
end
data_2D_b3 = [];
for d = 1:ndeltas
    for b=1:nshells
        for s=1:nshots
                data_2D_b3 = [data_2D_b3; squeeze(data_4D_b3(b, d, s, :)).'];
        end
    end
end
data_2D_b4 = [];
for d = 1:ndeltas
    for b=1:nshells
        for s=1:nshots
                data_2D_b4 = [data_2D_b4; squeeze(data_4D_b4(b, d, s, :)).'];
        end
    end
end
data_2D_b5 = [];
for d = 1:ndeltas
    for b=1:nshells
        for s=1:nshots
                data_2D_b5 = [data_2D_b5; squeeze(data_4D_b5(b, d, s, :)).'];
        end
    end
end

clear data d aux f fmax slice 

[nshells, ndeltas, nshots, npoints] = size(data_4D);

%% Quick check
[ft_single_4D, ft_sum_4D]                   = prepareft_4D(data_4D, nt, grpdly);
delta_ctr = 3;
figure
plot(ppmscale,real(squeeze(ft_sum_4D(1,delta_ctr,:))),'k')
hold on
plot(ppmscale,real(squeeze(ft_sum_4D(2,delta_ctr,:))),'r')
plot(ppmscale,real(squeeze(ft_sum_4D(3,delta_ctr,:))),'b')
plot(ppmscale,real(squeeze(ft_sum_4D(4,delta_ctr,:))),'m')
plot(ppmscale,real(squeeze(ft_sum_4D(5,delta_ctr,:))),'g')
clear delta_ctr

%% DENOISE

% Tensor MP-PCA 4D
[data_dn_tmp4D, sigma_tmp4D, P_tmp4D] = denoise_array_recursive_tensor(data_4D);

% Tensor MP-PCA 2D
[data_dn_tmp2D_aux, sigma_tmp2D, P_tmp2D] = denoise_array_recursive_tensor(data_2D);
data_dn_tmp2D = make_4D(data_dn_tmp2D_aux,data_4D,nt);

% Normal MP-PCA
[data_dn_mp_aux, sigma_mp, P_mp]    = denoiseMatrix(data_2D);
data_dn_mp = make_4D(data_dn_mp_aux,data_4D,nt);

% Patch Tensor MP-PCA 4D
data_dn_tmp4D_patch = denoise_patch_4D(data_4D,2,2,80,3000);

% Patch MP-PCA
data_dn_mp_patch_aux = denoise_patch_2D(data_2D,60, 3000,'MP-PCA');
data_dn_mp_patch = make_4D(data_dn_mp_patch_aux,data_4D,nt);

% Tensor MP-PCA 4D - one delta
[data_dn_tmp4D_DeltaMax, sigma_tmp4D_DeltaMax, P_tmp4D_DeltaMax] = denoise_array_recursive_tensor(data_4D_DeltaMax);
[data_dn_tmp4D_DeltaMin, sigma_tmp4D_DeltaMin, P_tmp4D_DeltaMin] = denoise_array_recursive_tensor(data_4D_DeltaMin);

% Tensor MP-PCA 4D - one bvalue
[data_dn_tmp4D_b1, sigma_tmp4D_b1, P_tmp4D_b1] = denoise_array_recursive_tensor(data_4D_b1);
[data_dn_tmp4D_b2, sigma_tmp4D_b2, P_tmp4D_b2] = denoise_array_recursive_tensor(data_4D_b2);
[data_dn_tmp4D_b3, sigma_tmp4D_b3, P_tmp4D_b3] = denoise_array_recursive_tensor(data_4D_b3);
[data_dn_tmp4D_b4, sigma_tmp4D_b4, P_tmp4D_b4] = denoise_array_recursive_tensor(data_4D_b4);
[data_dn_tmp4D_b5, sigma_tmp4D_b5, P_tmp4D_b5] = denoise_array_recursive_tensor(data_4D_b5);

data = [data_dn_tmp4D_b1;data_dn_tmp4D_b2;data_dn_tmp4D_b3;data_dn_tmp4D_b4;data_dn_tmp4D_b5];

save(fullfile(path_data,'sub-01_ses-01_combined_dmrs_allDelta-allb_FID_denoised.mat'),'data')

% MP-PCA - one bvalue
[data_dn_mp_aux, sigma_mp, P_mp]    = denoiseMatrix(data_2D_b5);
data_dn_mp_b5 = make_4D(data_dn_mp_aux,data_4D_b5,nt);
[data_dn_mp_aux, sigma_mp, P_mp]    = denoiseMatrix(data_2D_b1);
data_dn_mp_b1 = make_4D(data_dn_mp_aux,data_4D_b1,nt);

[nshells, ndeltas, nshots, npoints] = size(data_4D);

%% GET FID

% Get FID for normal denoising
[ft_single, ft_sum]                         = prepareft_4D(data_4D, nt, grpdly);
[ft_dn_tmp4D_single, ft_dn_tmp4D_sum]       = prepareft_4D(data_dn_tmp4D, nt, grpdly);
[ft_dn_tmp2D_single, ft_dn_tmp2D_sum]       = prepareft_4D(data_dn_tmp2D, nt, grpdly);
[ft_dn_mp_single, ft_dn_mp_sum]             = prepareft_4D(data_dn_mp, nt, grpdly);
[ft_dn_tmp4D_patch_single, ft_dn_tmp4D_patch_sum] = prepareft_4D(data_dn_tmp4D_patch, nt, grpdly);
[ft_dn_mp_patch_single, ft_dn_mp_patch_sum] = prepareft_4D(data_dn_mp_patch, nt, grpdly);


% Get FID for single delta denoising
[ft_single_DeltaMax, ft_sum_DeltaMax]           = prepareft_4D(data_4D_DeltaMax, nt, grpdly);
[ft_single_DeltaMin, ft_sum_DeltaMin]           = prepareft_4D(data_4D_DeltaMin, nt, grpdly);
[ft_dn_tmp4D_DeltaMax_single, ft_dn_tmp4D_DeltaMax_sum] = prepareft_4D(data_dn_tmp4D_DeltaMax, nt, grpdly);
[ft_dn_tmp4D_DeltaMin_single, ft_dn_tmp4D_DeltaMin_sum] = prepareft_4D(data_dn_tmp4D_DeltaMin, nt, grpdly);

% Get FID for single bvalue denoising
[ft_single_b5, ft_sum_b5]           = prepareft_4D(data_4D_b5, nt, grpdly);
[ft_single_b1, ft_sum_b1]           = prepareft_4D(data_4D_b1, nt, grpdly);
[ft_dn_tmp4D_b5_single, ft_dn_tmp4D_b5_sum] = prepareft_4D(data_dn_tmp4D_b5, nt, grpdly);
[ft_dn_tmp4D_b1_single, ft_dn_tmp4D_b1_sum] = prepareft_4D(data_dn_tmp4D_b1, nt, grpdly);
[ft_dn_mp_b5_single, ft_dn_mp_b5_sum] = prepareft_4D(data_dn_mp_b5, nt, grpdly);
[ft_dn_mp_b1_single, ft_dn_mp_b1_sum] = prepareft_4D(data_dn_mp_b1, nt, grpdly);


%% PLOT RESULTS
close all

% MPPCA vs tPMPCA 2D vs tMPPCA 4D
make_figure(ft_sum, ft_dn_tmp2D_sum,ft_dn_tmp4D_sum, ppmscale, bvals, deltas,...
    {'Original','t-MP-PCA 2D','res t-MP-PCA 2D','t-MP-PCA 4D','res t-MP-PCA 4D'},...
    fullfile(path_save,'Original_vs_tMPPCA2D_vs_tMPPCA4D_summed'))

make_figure(ft_sum, ft_dn_mp_sum,ft_dn_tmp4D_sum, ppmscale, bvals, deltas,...
    {'Original','MP-PCA','res MP-PCA','t-MP-PCA 4D','res t-MP-PCA 4D'},...
    fullfile(path_save,'Original_vs_MPPCA_vs_tMPPCA4D_summed'))

% Patch
make_figure(ft_sum, ft_dn_tmp4D_sum,ft_dn_tmp4D_patch_sum, ppmscale, bvals, deltas,...
    {'Original','t-MP-PCA 4D','res t-MP-PCA 4D','t-MP-PCA 4D patch','res t-MP-PCA 4D patch'},...
    fullfile(path_save,'Original_vs_tMPPCA4D_vs_tMPPCA4DPATCH_summed'))

make_figure(ft_sum, ft_dn_mp_patch_sum,ft_dn_tmp4D_patch_sum, ppmscale, bvals, deltas,...
    {'Original','MP-PCA patch','res MP-PCA patch','t-MP-PCA 4D patch','res t-MP-PCA 4D patch'},...
    fullfile(path_save,'Original_vs_MPPCAPATCH_vs_tMPPCA4DPATCH_summed'))

% Individual deltas
make_figure(ft_sum_DeltaMax, ft_dn_tmp4D_sum,ft_dn_tmp4D_DeltaMax_sum, ppmscale, bvals, deltas(end),...
    {'Original','t-MP-PCA 4D','res t-MP-PCA 4D','t-MP-PCA 4D 1Delta','res t-MP-PCA 4D 1Delta'},...
    fullfile(path_save,'Original_vs_together_vs_individualDeltaMax_summed'))
make_figure(ft_sum_DeltaMin, ft_dn_tmp4D_sum,ft_dn_tmp4D_DeltaMin_sum, ppmscale, bvals, deltas(1),...
    {'Original','t-MP-PCA 4D','res t-MP-PCA 4D','t-MP-PCA 4D 1Delta','res t-MP-PCA 4D 1Delta'},...
    fullfile(path_save,'Original_vs_together_vs_individualDeltaMin_summed'))

% Individual bvals
make_figure(ft_sum_b5, ft_dn_tmp4D_sum(end,:,:),ft_dn_tmp4D_b5_sum, ppmscale, bvals(end), deltas,...
    {'Original','t-MP-PCA 4D','res t-MP-PCA 4D','t-MP-PCA 4D b5','res t-MP-PCA 4D b5'},...
    fullfile(path_save,'Original_vs_together_vs_individual_tMPPCA_B5_summed'))
make_figure(ft_sum_b1, ft_dn_tmp4D_sum(1,:,:),ft_dn_tmp4D_b1_sum, ppmscale, bvals(1), deltas,...
    {'Original','t-MP-PCA 4D','res t-MP-PCA 4D','t-MP-PCA 4D b1','res t-MP-PCA 4D b1'},...
    fullfile(path_save,'Original_vs_together_vs_individual_tMPPCA_B1_summed'))

make_figure(ft_sum_b5, ft_dn_mp_sum(end,:,:),ft_dn_mp_b5_sum, ppmscale, bvals(end), deltas,...
    {'Original','MP-PCA','res MP-PCA','MP-PCA b5','res MP-PCA b5'},...
    fullfile(path_save,'Original_vs_together_vs_individual_MPPCA_B5_summed'))
make_figure(ft_sum_b1, ft_dn_mp_sum(1,:,:),ft_dn_mp_b1_sum, ppmscale, bvals(1), deltas,...
    {'Original','MP-PCA 4D','res MP-PCA 4D','MP-PCA 4D b1','resMP-PCA 4D b1'},...
    fullfile(path_save,'Original_vs_together_vs_individual_MPPCA_B1_summed'))


close all

%% GET RATIOS RESULTS

% MP-PCA vs t-MP-PCA 2D vs vs t-MP-PCA 4D
fprintf('********** \n MP-PCA vs t-MP-PCA 2D vs vs t-MP-PCA 4D\n')
compute_ratio_wdelta(ft_single, ppm_mask, 'raw data');
compute_ratio_wdelta(ft_dn_mp_single, ppm_mask, 'MP-PCA denoised data');
compute_ratio_wdelta(ft_dn_tmp2D_single, ppm_mask, 't-MP-PCA 2D denoised data');
compute_ratio_wdelta(ft_dn_tmp4D_single, ppm_mask, 't-MP-PCA 4D denoised data');

% Patch
fprintf('********** \n patch vs no patch \n')
compute_ratio_wdelta(ft_single, ppm_mask, 'raw data');
compute_ratio_wdelta(ft_dn_mp_single, ppm_mask, 'MP-PCA denoised data');
compute_ratio_wdelta(ft_dn_mp_patch_single, ppm_mask, 'MP-PCA patch denoised data');

compute_ratio_wdelta(ft_dn_tmp4D_single, ppm_mask, 't-MP-PCA 4D denoised data');
compute_ratio_wdelta(ft_dn_tmp4D_patch_single, ppm_mask, 't-MP-PCA 4D patch denoised data');

% Individual deltas
fprintf('********** \n individual vs together delta \n')
compute_ratio_wdelta(ft_single, ppm_mask, 'raw data');
compute_ratio_wdelta(ft_dn_tmp4D_single, ppm_mask, 't-MP-PCA 4D denoised data');
compute_ratio_wdelta(ft_dn_tmp4D_DeltaMax_single, ppm_mask, 't-MP-PCA 4D denoised data individual Delta Max');
compute_ratio_wdelta(ft_dn_tmp4D_DeltaMin_single, ppm_mask, 't-MP-PCA 4D denoised data individual Delta Min');

% Individual bvals
fprintf('********** \n individual vs together bval \n')
compute_ratio_wdelta(ft_single, ppm_mask, 'raw data');
compute_ratio_wdelta(ft_dn_tmp4D_single, ppm_mask, 't-MP-PCA 4D denoised data');
compute_ratio_wdelta_sb(ft_dn_tmp4D_b1_single,ft_dn_tmp4D_b5_single, ppm_mask, 't-MP-PCA 4D denoised data individual bval');
compute_ratio_wdelta_sb(ft_dn_mp_b1_single,ft_dn_mp_b5_single, ppm_mask, 'MP-PCA denoised data individual bval');


%% Auxiliar functions: prepareft_4D

function [ft_single, ft_sum] = prepareft_4D(data, nt, grpdly)

nshells  = size(data,1);
ndeltas  = size(data,2);
nshots   = size(data,3);

% SUM
for d = 1:ndeltas
    for s = 1:nshells
        fid = squeeze(data(s,d,:,:));
        fid=sum(fid(:,:),1);
        fid=[fid(round(grpdly(s))+1:end), zeros(1,round(grpdly(s)))];
        ft_sum(s,d,:)=fftshift(fft(fid./nt(s),[],2),2);
    end
end

% SINGLE
for d = 1:ndeltas
    for s = 1:nshells
        for shot=1:nshots

            fid = squeeze(data(s,d,shot,:))';
            fid=[fid(round(grpdly(s))+1:end), zeros(1,round(grpdly(s)))];
            ft_single(s,d,shot,:)=fftshift(fft(fid./nt(s),[],2),2);
        end
    end
end

end

%% Auxiliar functions: prepareft_2D
%
% function [ft_single, ft_sum] = prepareft_2D(data_2D, data_4D, nt, grpdly)
%
% nshells  = size(data_4D,1);
% ndeltas  = size(data_4D,2);
% nshots   = size(data_4D,3);
%
% idxs_shells = [0; cumsum(nt)];
% nd      = repmat(nshells, 1, ndeltas)';
% idxs_deltas = [0; cumsum(nd)].*nt;
%
% % SUM
% for d = 1:ndeltas
%     for s = 1:nshells
%
%         start_idx = idxs_deltas(d) + idxs_shells(s) + 1;
%         end_idx   = start_idx + nt(s)-1;
%
%         fid_block = squeeze(data_2D(start_idx:end_idx, :));
%
%         fid=sum(fid_block(:,:),1);
%         fid=[fid(round(grpdly(s))+1:end), zeros(1,round(grpdly(s)))];
%         ft_sum(d,s,:)=fftshift(fft(fid./nt(s),[],2),2);
%     end
% end
%
% % SINGLE
% for d = 1:ndeltas
%     for s = 1:nshells
%
%         start_idx = idxs_deltas(d) + idxs_shells(s) + 1;
%         end_idx   = start_idx + nt(s)-1;
%
%         for shot=1:nshots
%             fid = squeeze(data_2D(start_idx:end_idx, shot));
%             fid=[fid(round(grpdly(s))+1:end), zeros(1,round(grpdly(s)))];
%             ft_single(s,d,shot,:)=fftshift(fft(fid./nt(s),[],2),2);
%
%         end
%     end
%
% end
%
% end

%% Compute ratios

function ratio = compute_ratio_wdelta(ft, ppm_mask, mylabel)

ndeltas = size(ft,2);
nshots = size(ft,3);

fprintf('>> σ_bmin/σ_bmax in %s: \n', mylabel);

ratio_list = [];
for d = 1:ndeltas

    ratio=[];
    for shot=1:nshots

        ft_s_shot_min = squeeze(ft(1,d,shot,:));
        ft_s_shot_max = squeeze(ft(end,d,shot,:));
        sigma_bmin = std(real(ft_s_shot_min(ppm_mask)));
        sigma_bmax = std(real(ft_s_shot_max(ppm_mask)));
        ratio = [ratio; sigma_bmin/sigma_bmax];
    end

    fprintf('   for delta_indx=%d: %0.02f\n', d,mean(ratio));

    ratio_list = [ratio_list; mean(ratio)];

end

fprintf('   for mean delta: %0.02f\n', mean(ratio_list));
end

function ratio = compute_ratio_wdelta_sb(ft_b1,ft_b5, ppm_mask,mylabel )

fprintf('>> σ_bmin/σ_bmax in %s: \n', mylabel);
ratio_list = [];
for d = 1: size(ft_b1,2)

    ratio=[];
    for shot=1:size(ft_b1,3)

        ft_s_shot_min = squeeze(ft_b1(1,d,shot,:));
        ft_s_shot_max = squeeze(ft_b5(1,d,shot,:));
        sigma_bmin = std(real(ft_s_shot_min(ppm_mask)));
        sigma_bmax = std(real(ft_s_shot_max(ppm_mask)));
        ratio = [ratio; sigma_bmin/sigma_bmax];
    end

    fprintf('   for delta_indx=%d: %0.02f\n', d,mean(ratio));

    ratio_list = [ratio_list; mean(ratio)];

end

fprintf('   for mean delta: %0.02f\n', mean(ratio_list));
end
%% Make figure single
%
% function [] = make_figure_single(data1, data2, data3, ppmscale, legendlabel)
%
% % Get residuals
% for  shell=1:size(data1,2)
%     for shot=1:size(data1{shell},2)
%         data2_res{shell}{shot} = real(data2{shell}{shot})-real(data1{shell}{shot});
%         data3_res{shell}{shot} = real(data3{shell}{shot})-real(data1{shell}{shot});
%     end
% end
%
% % Define shot number to plot
% shot = 10;
%
% % Define number of shells
% nshells = size(data1,2);
%
% % Plot
% figure;
% tiledlayout(1, nshells, 'Padding', 'compact', 'TileSpacing', 'compact');
% set(gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
% for shell=1:nshells
%     subplot(1,nshells,shell)
%     plot(ppmscale,real(data1{shell}{shot}),'k')
%     hold on
%     plot(ppmscale,real(data2{shell}{shot}),'b')
%     plot(ppmscale,real(data2_res{shell}{shot})-700,'b--')
%     plot(ppmscale,real(data3{shell}{shot}),'r')
%     plot(ppmscale,real(data3_res{shell}{shot})-800,'r--')
%     title(sprintf('Shell %d, shot %d',shell,shot))
%     ylim([-1000 1000])
%     xlabel('ppm')
%     if shell==5
%         legend(legendlabel,'Location','south')
%     end
% end
%
% end

%% Make figure sum

function [] = make_figure(data1, data2, data3, ppmscale, bvals, deltas,legendlabel,outfile)

% Get residuals
data2_res = real(data2-data1);
data3_res = real(data3-data1);

% Define number of shells
nshells = size(data1,1);
ndeltas = size(data1,2);

% data2_res_std = std(reshape(data2_res,nshells, []), 0, 2);  
% data3_res_std = std(reshape(data3_res,nshells, []), 0, 2);  

% Plot
for d=1:ndeltas
    figure;
    tiledlayout(1, nshells, 'Padding', 'compact', 'TileSpacing', 'compact');
    set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 40, 20]);

    for shell=1:nshells
        nexttile
        plot(ppmscale,real(squeeze(data1(shell,d,:)))/1e6,'k')
        hold on
        plot(ppmscale,real(squeeze(data2(shell,d,:)))/1e6,'b')
        plot(ppmscale,real(squeeze(data2_res(shell,d,:)))/1e6-3,'b-')
        %
        plot(ppmscale,real(squeeze(data3(shell,d,:)))/1e6,'r')
        plot(ppmscale,real(squeeze(data3_res(shell,d,:)))/1e6-5,'r-')
        title(sprintf('b=%0.2f ms/{\\mu}m^2',bvals(shell)))
        ylim([-7 14])
        xlabel('ppm')
        ylabel('Signal×10^6')
        if shell==5
            legend(legendlabel,'Location','northeast')
        end
    end
    sgtitle(sprintf('Summed spectrum: \\Delta=%d ms',deltas(d)))
    saveas(gcf,sprintf('%s_Delta_%d_summed.png',outfile,deltas(d)))

end

end


%% Denoise patch_4D

function X_dn = denoise_patch_4D(X,patch_b,patch_d, patch_t,patch_p)

X_dn = zeros(size(X));               % Denoised output
counts = zeros(size(X));             % For overlap normalization

[nb, nd, nt, np] = size(X);
stride_b        = max([1, floor(patch_b / 3)]);
stride_d        = max([1, floor(patch_d / 3)]);
stride_np       = max([1, floor(patch_p / 3)]);
stride_nt       = max([1, floor(patch_t / 3)]);

for bi = 1:stride_b:nb
    for di = 1:stride_d:nd
        for ti = 1:stride_nt:nt
            for pi = 1:stride_np:np

                b_start = max(bi, 1);
                b_end   = min(bi + patch_b - 1, nb);

                d_start = max(di, 1);
                d_end   = min(di + patch_d - 1, nd);

                t_start = max(ti, 1);
                t_end   = min(ti + patch_t - 1, nt);

                p_start = max(pi, 1);
                p_end   = min(pi + patch_p - 1, np);

                patch = X(b_start:b_end, d_start:d_end, t_start:t_end, p_start:p_end);

                if ~anynan(patch)
                    % Run tensor MP-PCA
                    [patch_dn, ~, ~] = denoise_array_recursive_tensor(patch);

                    % Add to output and accumulate count for overlap
                    X_dn(b_start:b_end, d_start:d_end, t_start:t_end, p_start:p_end) = ...
                        X_dn(b_start:b_end, d_start:d_end, t_start:t_end, p_start:p_end) + patch_dn;

                    counts(b_start:b_end, d_start:d_end, t_start:t_end, p_start:p_end) = ...
                        counts(b_start:b_end, d_start:d_end, t_start:t_end, p_start:p_end) + 1;
                end
            end
        end
    end
end

% Normalize overlapping areas
X_dn = X_dn ./ max(counts, 1);
end


% %% Denoise patch_3D
% function X_dn = denoise_patch_3D(X,patch_b,patch_t,patch_f)
% 
% X_dn = zeros(size(X));               % Denoised output
% counts = zeros(size(X));             % For overlap normalization
% 
% [nb, nt, nf] = size(X);
% stride_b = max([1, floor(patch_b / 3)]);
% stride_t = max([1, floor(patch_t / 3)]);
% stride_f = max([1, floor(patch_f / 3)]);
% 
% for bi = 1:stride_b:nb
%     for ti = 1:stride_t:nt
%         for fi = 1:stride_f:nf
% 
%             b_start = max(bi, 1);
%             b_end   = min(bi + patch_b - 1, nb);
% 
%             t_start = max(ti, 1);
%             t_end   = min(ti + patch_t - 1, nt);
% 
%             f_start = max(fi, 1);
%             f_end   = min(fi + patch_f - 1, nf);
% 
%             patch = X(b_start:b_end, t_start:t_end, f_start:f_end);
% 
%             if ~anynan(patch)
%                 % Run tensor MP-PCA
%                 [patch_dn, ~, ~] = denoise_array_recursive_tensor(patch);
% 
%                 % Add to output and accumulate count for overlap
%                 X_dn(b_start:b_end, t_start:t_end, f_start:f_end) = ...
%                     X_dn(b_start:b_end, t_start:t_end, f_start:f_end) + patch_dn;
% 
%                 counts(b_start:b_end, t_start:t_end, f_start:f_end) = ...
%                     counts(b_start:b_end, t_start:t_end, f_start:f_end) + 1;
%             end
%         end
%     end
% end
% 
% % Normalize overlapping areas
% X_dn = X_dn ./ max(counts, 1);
% end

function data_4D = make_4D(data_2D,data_4D_dim,nt)
    nshells  = size(data_4D_dim,1);
    ndeltas  = size(data_4D_dim,2);
    nshots   = size(data_4D_dim,3);
    npoints   = size(data_4D_dim,4);

    idxs_shells = [0; cumsum(nt)];
    nd      = repmat(nshells, 1, ndeltas)';
    idxs_deltas = [0; cumsum(nd)]*nt(1);
    data_4D=zeros(nshells,ndeltas,nshots,npoints);
    for d = 1:ndeltas
        for b=1:nshells
            start_idx = idxs_deltas(d) + idxs_shells(b) + 1;
            end_idx   = start_idx + nt(b)-1;
            data_4D(b,d,:,:) = squeeze(data_2D(start_idx:end_idx, :));
        end
    end
end


%% Denoise patch_2D
function X_dn = denoise_patch_2D(X, patch_bt,patch_f, type)
%patch_bt = 50;      % Number of b-values and shots per patch (e.g., 3)
%patch_f = 1000;    % Number of FID points per patch
X_dn = zeros(size(X));               % Denoised output
counts = zeros(size(X));             % For overlap normalization

[nb, nf] = size(X);
stride_bt        = max([1, floor(patch_bt / 2)]);
stride_f           = max([1, floor(patch_f / 2)]);

X(X==0)=NaN;

for bti = 1:stride_bt:nb
    for fi = 1:stride_f:nf

        bt_start = max(bti, 1);
        bt_end   = min(bti + patch_bt - 1, nb);

        f_start = max(fi, 1);
        f_end   = min(fi + patch_f - 1, nf);

        patch = X(bt_start:bt_end,f_start:f_end);

        if ~anynan(patch)
            % Run tensor MP-PCA (you already have this)
            if strcmp(type,'MP-PCA')
                [patch_dn, ~, ~] = denoiseMatrix(patch);
            elseif strcmp(type,'t-MP-PCA')
                [patch_dn, ~, ~] = denoise_array_recursive_tensor(patch);
            end
            % Add to output and accumulate count for overlap
            X_dn(bt_start:bt_end, f_start:f_end) = ...
                X_dn(bt_start:bt_end, f_start:f_end) + patch_dn;
    
            counts(bt_start:bt_end, f_start:f_end) = ...
                counts(bt_start:bt_end, f_start:f_end) + 1;
        end
    end
end


% Normalize overlapping areas
X_dn = X_dn ./ max(counts, 1);
end
