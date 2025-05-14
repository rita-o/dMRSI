function denoise_in_matlab(inputs,output,header_f,N)
% Function that uses tensor MP-PCA implementation of matlab to denoise the
% dwi data containing multiple diffusion times
% Data = [x, y, z, bvals/bvecs, diffTime]
% To ba called from python script
%
% INPUTS: 
%       inputs - cell array {1xM} with path to the files containing dwi 4D
%           datasets of a give diffusion time, it has M diffusion times
%       output - string to path where to save the denoised data. The data
%           will be concatenated across diffusion times, similarly to other
%           methods and for the easy of use afterwards
%       header_f - string to path where the concatenated data (across diffusion 
%           times is). Serves as header to save the denoised data.
%       N  - patch/kernel size to use when denoising
% OUTPUTS:
%       none
% 
% __________________________
% Rita Oliveira, 
% Mai 2025,
% MicMap Lab, Switzerland

    addpath(genpath('/home/localadmin/Documents/Rita/Toolboxes/spm12')) 
    addpath(genpath('/home/localadmin/Documents/Rita/Toolboxes/Tensor-MP-PCA-main'))
    disp('#########################')
    disp('>>> We are in matlab now')

    % Build 5D matrix with diffusion times in the 5th dimension
    for i=1:length(inputs)
        data_5D(:,:,:,:,i) = spm_read_vols(spm_vol(inputs{i}));
    end

    % % Convert input string to number
    N = str2num(N);

    % Denoised 5D matrix
    disp('#########################')
    disp('>>> Denoising data, it might take a while...')
    [denoised_5D, sigma_5D, P_5D] = denoise_recursive_tensor(data_5D,[N N N]);

    % Save denoised data
    disp('#########################')
    disp('>>> Saving denoised data')
    header = spm_vol(header_f);
    s = size(denoised_5D);  % s = [80, 22, 64, 126, 3]
    for i=1:s(4)*s(5)
         header(i).fname  = output;
         fprintf('saving volume %d \n',i)
         spm_write_vol(header(i),squeeze(denoised_5D(:,:,:,i)));
    end
  

end
