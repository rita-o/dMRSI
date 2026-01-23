function processing_dmrs_matlab(input_path, output_path, data_folders_list, coil_type, toolbox_path, basis_set, LCMpath)

% Add paths
addpath(genpath(toolbox_path)) 

% Define folders
foldersum     = fullfile(output_path, 'processed','sum');
folderrawsave = fullfile(output_path, 'raw');
folder_quanti_save = fullfile(output_path,'quantified');

% Loop through data folders
data_folders_list = str2num(data_folders_list);

for j = 1:length(data_folders_list)

        % Get scan id
        scan_id = data_folders_list(j);
        [~, raw_data,~] = fileparts(input_path);
        disp('****************************************************************')
        fprintf('Analysing data \n %s \n with scan number %d ... \n \n',raw_data, scan_id)

        % Get data from Bruker
        fprintf('   Get data from Bruker by execution a_Create_study_Bruker_9p4T_EM (Toi version 10.11.2025). \n')
        a_Create_study_Bruker_RatCryo_9p4T_EM_MB(input_path, output_path, scan_id, coil_type); % rat cryo

        % Process data
        fprintf('  Process data by execution b_Processing_FidA_9p4T_EM (Toi version 10.11.2025). \n')
        phi = 5.6;% [0; 6.3]
        b_Preprocessing_FidA_9p4T_EM_MB(output_path , scan_id, phi);
        fprintf('\n Finished! \n')
        disp('****************************************************************')
        close all;
    
end

% Quantitfy
disp('****************************************************************')
fprintf('Quantify spectrum with LC Model for all data folders ... \n')
Quantif_LCModel_noabsquant_MB(foldersum, folderrawsave, folder_quanti_save,basis_set,LCMpath, output_path)

end