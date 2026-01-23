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
        
        fprintf('\n');
        disp('==============================================================');
        fprintf(' Processing dataset: %s\n', raw_data);
        fprintf(' Scan number       : %d\n', scan_id);
        disp('==============================================================');

        % ------------------------------------------------------------------
        % Import data from Bruker
        fprintf('\n[1/2] Importing Bruker data\n');
        fprintf('      with a_Create_study_Bruker_RatCryo_9p4T_EM_MB\n');
        fprintf('        (Toi version 10.11.2025)\n\n');
    
        a_Create_study_Bruker_RatCryo_9p4T_EM_MB( ...
            input_path, output_path, scan_id, coil_type);

        % ------------------------------------------------------------------
        % Preprocessing
        fprintf('\n[2/2] Preprocessing data\n');
        fprintf('      with b_Preprocessing_FidA_9p4T_EM_MB\n');
        fprintf('        (Toi version 10.11.2025)\n\n');
    
        phi = 5.6;  % phase correction
        b_Preprocessing_FidA_9p4T_EM_MB(output_path, scan_id, phi);
    
        fprintf('\n ✔ Finished scan %d\n', scan_id);
        disp('--------------------------------------------------------------');
    
end

% ----------------------------------------------------------------------
% Quantification
fprintf('\n');
disp('==============================================================');
fprintf(' Quantifying spectra with LCModel for all datasets\n');
disp('==============================================================');

Quantif_LCModel_noabsquant_MB( ...
    foldersum, folderrawsave, folder_quanti_save, ...
    basis_set, LCMpath, output_path);

end