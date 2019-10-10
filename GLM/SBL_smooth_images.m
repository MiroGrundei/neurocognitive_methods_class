%% SBL: smooth nifti scalptime image
% -------------------------------------------------------------------------
clc; close; clear all;
addpath('C:\Users\Miro\Documents\Code\MATLAB\toolboxes\spm12')
spm('defaults', 'EEG');
spm_jobman('initcfg');
project_dir = fullfile('D:','EEG','SCAN_class');

%% 

smoothing_kernel = [12 12 0];
prefix = 'image_bletdhmar_';
img_name = 'condition_1.nii';
SJs = [32,33];

% loop through subjects
for s = SJs 
    
    subject = sprintf('sub-%02d',s);
    fprintf('Smoothing image of %s \n',subject)                             % user info
    
    fname = fullfile(project_dir,'data',subject,'preprocessed',             ...
                     [prefix subject '_SBL'], img_name);                    % converted nifti file
    
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(fname);        
    matlabbatch{1}.spm.spatial.smooth.fwhm = smoothing_kernel;
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 1;
    matlabbatch{1}.spm.spatial.smooth.prefix = sprintf('s%dx%dx%d_',        ...
                                                       smoothing_kernel(1), ...
                                                       smoothing_kernel(2), ...
                                                       smoothing_kernel(3));

    spm_jobman('run', matlabbatch);                                         % run
    clear matlabbatch
    
end

disp('Done!')