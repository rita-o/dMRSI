%script_for_Pauline
%
%  !!!!!!!!!!!! EXECUTER LE SCRIPT DANS LE BON DOSSIER !!!!!!!!!!!!!!!!
%

%% GLOBAL FOLDERS
addpath(genpath('/data/u_mbrammerloh_diss/new_owncloud/projects/CHUV/Data_Pauline/Code_quickProcessing_dMRS_Pauline/dSPECIAL_Manual_Codes_Toi')) % RITA: update this path
% addpath('C:\Users\mouge/home/leal/Desktop/cudalbu_biblioFP/Code_matlabl\Documents\MATLAB\codes\codes\Codes_CIBM\dSPECIAL_Manual_Codes_Toi\FID-A-master')

%%  STEP 0: (OPTIONAL) CLEAR/CLOSE EVERYTHING
clear all 
close all 
clc

%% STEP 1: CHANGE PARAMETERS FOR EXPERIMENT
% (CTRL+ ENTREE)
data_folder = '/data/u_mbrammerloh_diss/new_owncloud/projects/CHUV/Data_Pauline';% RITA: update this path

rat_number = 1778;
raw_folder_name = '20251114_090217_20251111_CTD_P10_PL_20251114_Rat_1778_M_CTD_1_6';

%coil_type = 'rat'; % for P30, rats scanned with rat cryo probe
coil_type = 'mouse'; % for P10, rat pups scanned with moise cryo probe

expnb_min = 17; % SCAN ID
expnb_max = 46; % 54;
expnb_to_skip = [24 25 33 34 35 36 42];

expnb_switch_basis_set = [35 46]; % this is still ignored :/

% derived folder names
folder_exp =string([data_folder '/' num2str(rat_number) '/raw_data/' raw_folder_name '/']); % PATH SCAN
folder_results = ['/data/u_mbrammerloh_diss/new_owncloud/projects/CHUV/Data_Pauline/' num2str(rat_number) '/derived/mrs4brain/']; % RITA: update this path

basis_sets = { ...
    '/data/u_mbrammerloh_diss/new_owncloud/projects/CHUV/dmrs_basis/Basis_Set_dSPECIAL_differentTM/9p4T_Toi_dSPECIAL_TM250_20260106.BASIS', ...
    '/data/u_mbrammerloh_diss/new_owncloud/projects/CHUV/dmrs_basis/Basis_Set_dSPECIAL_differentTM/9p4T_Toi_dSPECIAL_TM150_20260106.BASIS', ...
    '/data/u_mbrammerloh_diss/new_owncloud/projects/CHUV/dmrs_basis/Basis_Set_dSPECIAL_differentTM/9p4T_Toi_dSPECIAL_TM41_20250925.BASIS' ...
    };% RITA: update this path

LCMpath='/data/u_mbrammerloh_diss/tools/LCModel/binaries/linux/'; % RITA: update this path
foldersum=[folder_results, 'processed/sum/'];
folderrawsave=[folder_results,'raw/'];
folder_quanti_save = [folder_results,'quantified/'];

for expnb = expnb_min:expnb_max
    if ~ismember(expnb, expnb_to_skip)
        basis_set = basis_sets{1};
        for i = 1:length(expnb_switch_basis_set)
            if expnb >= expnb_switch_basis_set(i)
                basis_set = basis_sets{i+1};
            end
        end

        %a_Create_study_Bruker_MouseCryo_9p4T_EM(folderexp , expnb); % mouse cryo
        a_Create_study_Bruker_RatCryo_9p4T_EM_MB(folder_exp, folder_results, expnb, coil_type); % rat cryo
        phi = 5.6;% [0; 6.3]
        b_Preprocessing_FidA_9p4T_EM_MB(folder_results , expnb, phi);
        close all;
    end
end

basis_set = basis_sets{1};
Quantif_LCModel_noabsquant_MB(foldersum, folderrawsave, folder_quanti_save,basis_set,LCMpath, folder_results)
