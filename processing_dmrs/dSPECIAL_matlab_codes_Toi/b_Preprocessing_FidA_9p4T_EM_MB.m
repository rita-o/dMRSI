
function [] = b_Preprocessing_FidA_9p4T_EM_MB(folder_results , expnb, phi)
%% b_Preprocessing_FidA_9p4T_EM
%EM no change from Toi except path version 10.11.2025
% RO Rita changed absolute paths to relative paths and used fullfile
% instead of hardcoded slashes
% clear; clc; %close all; 
% 
% folder_results   = "Z:\DATA\CSI\9.4T-FIDCSI-Data\20250424_084847_Toi_dMRSI_20250424_RatFem_dMRI_dSPECIAL_dMRSI_1_11_cc";
% 
% expnb       = 35; % scan number
%disp(['execution b_Processing_FidA_9p4T_EM (Toi version 10.11.2025) on data: ',folder_results ])
tic
%% step 1 - do preprocessing with FIDA
filelist = dir(fullfile(folder_results, "raw","*ser.mat"));

idx = find(contains({filelist.name}, ['_',num2str(expnb),'_']));   % look for number 43 in name
filelist = filelist(idx,:); 

% add FID-A master to Matlab path
% addpath(genpath('Y:\Toi\dMRS_SPECIAL\dSPECIAL_scripts\dSPECIAL_Manual_Codes\FID-A-master'))

LBall               = 10;
doprocessing        = true;
dosavesum           = true;
dosaveprocessing    = true; 

spectralreginf      = 0.5; %ppm
spectralregsup      = 4.3; %ppm

if doprocessing
screenSz = get(0,'ScreenSize'); figure('Position',[0,screenSz(4)/3-150,screenSz(3)-screenSz(3)/5,screenSz(4)/3],'Color','w');
til = tiledlayout(length(filelist),4,"Padding","tight");

for file=1:length(filelist)
    disp(filelist(file).name)
    % f1=figure;

    %% 1-Convert Bruker study to FID A structure
    out0a=convertNicoStudyFidAformat_isison(char(fullfile(folder_results,'raw')), filelist(file).name);
    out0b=convertNicoStudyFidAformat_isisoff(char(fullfile(folder_results,'raw')), filelist(file).name);
    %apply LB
    dw=out0a.dwelltime;
    tt=[0:dw:dw*(out0a.n-1)];
    out0alb=out0a;
    fids0alb=out0alb.fids.*repmat(exp(-tt*pi*LBall).',1,out0a.averages);
    out0alb.fids=fids0alb;
    out0alb.specs=fftshift(fft(out0alb.fids.',[],2),2).';
    out0blb=out0b;
    fids0blb=out0blb.fids.*repmat(exp(-tt*pi*LBall).',1,out0b.averages);
    out0blb.fids=fids0blb;
    out0blb.specs=fftshift(fft(out0blb.fids.',[],2),2).';
    

    % figure(f1);
    % subplot(2,2,1)
    nexttile(til);
    for k=1:out0a.averages
        plot(real(out0a.specs(:,k)))
        hold on; xlim tight;  ylim  tight;
        % xlim([2900,3050])
    end
    for k=1:out0b.averages
        plot(real(out0b.specs(:,k)))
        hold on; xlim tight;  ylim  tight;
        % xlim([2900,3050])
    end 
    title(['E',num2str(expnb),' - raw'])
    
    % figure(f1);
    % subplot(2,2,2)
    nexttile(til);
    for k=1:out0alb.averages
        plot(real(out0alb.specs(:,k)))
        hold on; xlim tight;  ylim  tight;
        % xlim([2900,3050])
    end 
    for k=1:out0blb.averages
        plot(real(out0blb.specs(:,k)))
        hold on; xlim tight;  ylim  tight;
        % xlim([2900,3050])
    end 
    title(['E',num2str(expnb),' - Linebroadening ' num2str(LBall) 'Hz'])
    
    %% 2-align av.
    [out1alb,fsa,phsa]=op_alignAverages_fd_jm_conj(out0alb,4.7+spectralreginf,4.7+spectralregsup,0.5,'y');%5.2,9,0.5,'y');
    [out1blb,fsb,phsb]=op_alignAverages_fd_jm_conj(out0blb,4.7+spectralreginf,4.7+spectralregsup,0.5,'y');%5.2,9,0.5,'y');
    % fs        = Vector of frequency shifts (in Hz) used for alignment.
    % phs       = Vector of phase shifts (in degrees) used for alignment.

    % figure(f1);
    % subplot(2,2,3)
    nexttile(til);
    for k=1:out1alb.averages
        plot(real(out1alb.specs(:,k)))
        hold on; xlim tight; ylim  tight;
        % xlim([2900,3050])
    end
    for k=1:out1blb.averages
        plot(real(out1blb.specs(:,k)))
        hold on; xlim tight; ylim  tight;
        % xlim([2900,3050])
    end
    title(['E',num2str(expnb),' - Spectral alignment (JNear averaged)'])
    
    %remove LB
    out1a=out1alb;
    fids1a=out1alb.fids.*repmat(exp(tt*pi*LBall).',1,out1alb.averages);
    out1a.fids=fids1a;
    out1a.specs=fftshift(fft(out1a.fids.',[],2),2).';
    out1b=out1blb;
    fids1b=out1blb.fids.*repmat(exp(tt*pi*LBall).',1,out1blb.averages);
    out1b.fids=fids1b;
    out1b.specs=fftshift(fft(out1b.fids.',[],2),2).'; 

    %% 3-outlier removal 
    [out2a,metrica,badAveragesa]=op_rmbadaverages_jm(out1a,1.5,'f'); %performs 10Hz LB inside
    [out2b,metricb,badAveragesb]=op_rmbadaverages_jm(out1b,1.5,'f'); %performs 10Hz LB inside
    
    %apply LB
    out2alb=out2a;
    fids2alb=out2a.fids.*repmat(exp(-tt*pi*LBall).',1,out2a.averages);
%     out2alb.fids=conj(fids2alb);
    out2alb.fids=fids2alb;
    out2alb.specs=fftshift(fft(out2alb.fids.',[],2),2).';
    out2blb=out2b;
    fids2blb=out2b.fids.*repmat(exp(-tt*pi*LBall).',1,out2b.averages);
    out2blb.fids=fids2blb;
    out2blb.specs=fftshift(fft(out2blb.fids.',[],2),2).';
    
    % figure(f1);
    % subplot(2,2,4)
    nexttile(til);
    for k=1:out2alb.averages
        plot(real(out2alb.specs(:,k)))
        hold on; xlim tight;  ylim  tight;
        % xlim([2900,3050])
    end
    for k=1:out2blb.averages
        plot(real(out2blb.specs(:,k)))
        hold on; xlim tight;  ylim  tight;
        % xlim([2900,3050])
    end
    title(['E',num2str(expnb),' - Outlier removal']); ylim  tight;
    a=legend(gca,num2str([badAveragesa;badAveragesb]),'Box','off','Location','southwest'); title(a,'Bad averages');

    %combine on/off
    clear fidtot;
    fida=out1a.fids.'; 
    fidb=out1b.fids.'; 
    fidtot(1:2:size(fida,1)*2,:)=fida; 
    fidtot(2:2:size(fida,1)*2,:)=fidb;
    
    ind=1;
    clear fid2sum;
    for k=1:size(fidtot,1)/2
        fid2sum(ind,:)=sum(fidtot((k-1)*2+1:k*2,:)); 
        ind=ind+1;
    end 
    
    ind=1;
    clear fidmocor;
    for k=1:size(fidtot,1)/2 %from 1 to 80 pairs
        if ismember(k,badAveragesa) % if in pair k, the odd is an outlier (= present in the 1st outlier list) - out
        else % if in pair k, the odd is an not an outlier 
            if ismember(k,badAveragesb) % if in pair k, the even is an outlier (= present in the 2nd outlier list) - out
            else % if in pair k, neither the odd nor the even are outliers
                fidmocor(ind,:) = fid2sum(k,:);
                ind=ind+1; 
            end 
        end 
    end 

    fidmocor=conj(fidmocor);
    
    %% 4-add all the info to the Matlab study structure
    if dosaveprocessing
        load(char(fullfile(folder_results,'raw',filelist(file).name)));
        study.fidaprocess.phsa=phsa;
        study.fidaprocess.fsa=fsa;
        study.fidaprocess.metrica=metrica;
        study.fidaprocess.badAveragesa=badAveragesa; 
        study.fidaprocess.phsb=phsb;
        study.fidaprocess.fsb=fsb;
        study.fidaprocess.metricb=metricb;
        study.fidaprocess.badAveragesb=badAveragesb;

        study.params.nt=size(fidmocor,1)*2; 
        study.multiplicity=size(fidmocor,1);

        study.process.apodparam1=zeros(1,size(fidmocor,1));
        study.process.apodparam2=zeros(1,size(fidmocor,1));
        study.process.phasecorr0=zeros(1,size(fidmocor,1));
        study.process.phasecorr1=zeros(1,size(fidmocor,1));

        study.data.real=zeros(size(fidmocor,1),1,size(fidmocor,2));
        study.data.real(:,1,:)=real(fidmocor);
        study.data.imag=zeros(size(fidmocor,1),1,size(fidmocor,2));
        study.data.imag(:,1,:)=imag(fidmocor);

        if ~exist(fullfile(folder_results, 'processed'), 'dir')
           mkdir(fullfile(folder_results, 'processed'))
        end
        save(fullfile(folder_results, 'processed', [filelist(file).name(1:end-4), '_processed.mat']),'study')   
    end 
    end 
end 
    
%% STEP 2 - APPLY PREPROCESSING, SUM AND PHASE THE SUM 

if dosavesum

screenSz = get(0,'ScreenSize'); figure('Position',[screenSz(3)/5*4,screenSz(4)/3-150,screenSz(3)/5,screenSz(4)/3],'Color','w');

for file=1:length(filelist)
    load(fullfile(folder_results, 'processed', [filelist(file).name(1:end-4), '_processed.mat']),'study')   

    fidmocor=squeeze(study.data.real)+1i*squeeze(study.data.imag);
    sumfid=sum(fidmocor);%./(size(fidmocor,1).*2); 
   
    
    %% save
    study.data.real=zeros(1,1,study.np/2);
    study.data.imag=zeros(1,1,study.np/2);

    study.data.real(1,1,:)=real(sumfid); 
    study.data.imag(1,1,:)=imag(sumfid);

    study.multiplicity=1;
    study.process.lsfid=0;
    study.process.apodparam1=0;
    study.process.apodparam2=0;
    study.process.phasecorr0=0;
    study.process.phasecorr1=0;
    study.process.B0=zeros(1,study.np/2);

    filename=filelist(file).name(1:end-4);
    if ~exist(fullfile(folder_results, 'processed','sum'), 'dir')
        mkdir(fullfile(folder_results, 'processed','sum'))
    end
    save(fullfile(folder_results, 'processed','sum',['SUM_' filename '_processed.mat']),'study');
    
    %plot

    % phi = 0.7; %%PHASE 
        receiveroffset_ppm = 4.7;
        sfrq1H = study.resfreq;
        receiveroffset_hz=receiveroffset_ppm*sfrq1H;
        frqEch = study.spectralwidth;
        np = study.np/2;
        dw=1/frqEch;
        sifactor=1;
        fmax=1/(2*dw);
        f_vec=[-fmax:2*fmax/(sifactor*np-1):fmax];
        ppm_vec=(-f_vec+receiveroffset_hz)/sfrq1H;
        % time=((0:np-1)*dw)';


    % figure; 
    ftcorr=fftshift(fft(sumfid.*exp(-1i*phi)./study.params.nt,[],2),2); 
    plot(ppm_vec,real(ftcorr))
    hold on 
    load(char(fullfile(folder_results,'raw',[filelist(file).name(1:end-4) '.mat'])),'study')   
    fidini=squeeze(study.data.real)+1i*squeeze(study.data.imag);
    sumfidini=sum(fidini)./study.params.nt;
    sumfidini=[sumfidini(77:end),zeros(1,76)];  
    ftini=fftshift(fft(sumfidini.*exp(-1i*phi),[],2),2); 
    plot(ppm_vec,real(ftini)); xlim(gca,[2000 3500]);
    set(gca,'XDir','reverse')
    xlabel('\delta [ppm]')
    legend('Corrected','Original','Box','off','Location','southeast'); xlim tight;  ylim  tight;

end
end 

%disp('**************** done ****************')
toc
end