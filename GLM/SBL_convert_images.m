%% SBL: convert meeg file to nifti 
% -------------------------------------------------------------------------
clc; close; clear all;
addpath('C:\Users\Miro\Documents\Code\MATLAB\toolboxes\spm12')
spm('defaults', 'EEG');
spm_jobman('initcfg');
project_dir = fullfile('D:','EEG','SCAN_class');

%% Create images for each trial for each subject from preprocessed meeg file
% -------------------------------------------------------------------------
prefix = 'bletdhmar_';
SJs = [42];

% loop through subjects
for s = SJs 
    
    subject = sprintf('sub-%02d',s);
    fprintf('Creating .nii for %s \n',subject)                              % user info
    
    fname = fullfile(project_dir,'data',subject,'preprocessed',             ...
                     [prefix subject '_SBL.mat']);                          % preprocessed meeg file
    
    % remove unwanted trials before image conversion             
                 
    matlabbatch{1}.spm.meeg.images.convert2images.D = cellstr(fname);
    matlabbatch{1}.spm.meeg.images.convert2images.mode = 'scalp x time';    % spatiotemporal data
    matlabbatch{1}.spm.meeg.images.convert2images.conditions = {};          % use all conditions 
    matlabbatch{1}.spm.meeg.images.convert2images.channels{1}.type = 'EEG'; % channels to be used
    matlabbatch{1}.spm.meeg.images.convert2images.timewin = [-Inf Inf];     % if applicable: restrict time window of interest
    matlabbatch{1}.spm.meeg.images.convert2images.freqwin = [-Inf Inf];     % frequency window of interest
    matlabbatch{1}.spm.meeg.images.convert2images.prefix = 'image_';        % prefix for newly created .nii file
    
    spm_jobman('serial',matlabbatch);                                       % run the job
    clear matlabbatch
    
end

disp('Done!')
