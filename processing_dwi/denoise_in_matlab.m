function denoise_in_matlab(input,output,numvols,N, toolbox_path)
% Function that uses tensor MP-PCA implementation of matlab to denoise the
% dwi data containing multiple diffusion times. It will denoise a matrix in
% the form:
% Data = [x, y, z, bvals/bvecs, diffTime]
% To be called from python script.
%
% INPUTS: 
%       input - string to path where the concatenated data (across diffusion times) is.
%       output - string to path where to save the denoised data. 
%       numvols - vector (1xM) with number of volumes for each diffusion
%       time M. Used to split the input 4D volume (x,y,z,M*bvals/bvecs) 
%       into a 5D matrix (x,y,z,bvals/bvecs,M)
%       N  - patch/kernel size to use when denoising
%       toolbox_path - path to the folder where Matlab toolboxes are
%
% OUTPUTS:
%       none
% 
% __________________________
% Rita Oliveira, 
% Mai 2025,
% MicMap Lab, Switzerland

    addpath(genpath(fullfile(toolbox_path,'spm12'))) 
    addpath(genpath(fullfile(toolbox_path,'Tensor-MP-PCA-main'))) 
    disp('#########################')
    disp('>>> We are in matlab now')

    % Convert input strings to numbers
    N = str2num(N);
    numvols = str2num(numvols);

    % Build 5D matrix with diffusion times in the 5th dimension
    dwi_total = spm_read_vols(spm_vol(input));
    for i=1:length(numvols)
        data_5D(:,:,:,:,i) = dwi_total(:,:,:,numvols(i)*(i-1)+1:numvols(i)*i);
    end

    % Denoised 5D matrix
    disp('#########################')
    disp('>>> Denoising data, it might take a while...')
    [denoised_5D, sigma_5D, P_5D] = denoise_recursive_tensor(data_5D,[N N N]);

    % Save denoised data
    disp('#########################')
    disp('>>> Saving denoised data')
    header = spm_vol(input);
    s = size(denoised_5D);  
    for i=1:s(4)*s(5)
         header(i).fname  = output;
         fprintf('saving volume %d \n',i)
         spm_write_vol(header(i),squeeze(denoised_5D(:,:,:,i)));
    end
  

end
