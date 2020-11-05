%% Distributed (3D) source reconstruction with SPM

clc; clear; close

%% Set Paths
% -------------------------------------------------------------------------

addpath(genpath('C:\Users\Miro\Documents\Code\MATLAB\toolboxes\spm12'))
% spm('defaults', 'EEG');

base_dir = 'E:\SCAN_class\neurocognitive_methods_class\';                   % change to directory          
prefix = 'SEP_rbletdhmar_';                                                 % prefix from averaged and preprocessed meeg file
out_fol = 'distributed_SR';

subjects = [3,4];
% time windows of interest
wois = [18 23; ...                                                          % N20
        40 60; ...                                                          % P50
        110 160];                                                           % N140
    
%% Forward model inversion
% -------------------------------------------------------------------------

% loop through subjects
for s = subjects

    sub = sprintf('sub-%02d',s);
    fprintf('Processing subject: %s\n',sub)
    
    % dircetory for saving resulting file
    out_dir = fullfile(base_dir, 'data', 'source_recon', sub, out_fol);
    if ~exist(out_dir,'dir')
        mkdir(out_dir) 
    end 
    
    % load data & copy to out directory (dublicate)
    data = fullfile(base_dir,'data','erp',sub,sprintf('%s%s_SBL.mat',prefix,sub));
    D = spm_eeg_load(data);    
    S = [];
    S.D = D;
    S.outfile = fullfile(out_dir,fname(D));
    D = spm_eeg_copy(S);
    data = fullfile(out_dir,fname(D));
    
    % remove any existing inversions
    if isfield(D,'inv')
        D = rmfield(D,'inv');
    end

    % specify forward model
    clear matlabbatch
    matlabbatch{1}.spm.meeg.source.headmodel.D = cellstr(data);
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;                       % index of the model (if more than one saved in a file)
    matlabbatch{1}.spm.meeg.source.headmodel.comment = 'SEP';           
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;   % use headshape template
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;           % resolution of the mesh (of head template)
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas'; % use fiducials specified within the ERP file 
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';       % Boundary Element Method for head model
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    
    % run the job: create forward model
    spm_jobman('run', matlabbatch)

    % loop through time windows of interest
    for t = 1:size(wois,1)

        clear matlabbatch
        matlabbatch{1}.spm.meeg.source.invert.D = cellstr(data);
        matlabbatch{1}.spm.meeg.source.invert.val = t;                      % index of the inversion (here more than one saved in one file)
        matlabbatch{1}.spm.meeg.source.invert.whatconditions.condlabel = {'SEP'};   % condition to be used from ERP file
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS';
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.woi = wois(t,:); % time window of interest
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.foi = [0 256]; % if applicable: frequency of interest
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.hanning = 1; % hanning window
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''}; % spatial mask for prior locations
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3); % restrict source by prior locations Xby3 location array
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''}; % spatial mask to restrict analysis
        matlabbatch{1}.spm.meeg.source.invert.modality = {'EEG'};
        
        % run the job: invert forward model
        spm_jobman('run', matlabbatch)
    end

end

%% Create Nifti from inversion results (to work with SPM standard pipeline)
% -------------------------------------------------------------------------

% loop through subjects
for s = 1:numel(subjects)  
    
    sub = sprintf('sub-%02d',subjects(s));
    
    % get the data with inverted forward model
    data = fullfile(base_dir, 'data', 'source_recon', sub, out_fol,sprintf('%s%s_SBL.mat',prefix,sub));
    
    % loop through time windows of interest
    for t = 1:size(wois,1) 

        clear matlabbatch
        matlabbatch{1}.spm.meeg.source.results.D = cellstr(data);
        matlabbatch{1}.spm.meeg.source.results.val = t;
        matlabbatch{1}.spm.meeg.source.results.woi = [wois(t,:)];
        matlabbatch{1}.spm.meeg.source.results.foi = [0 0];
        matlabbatch{1}.spm.meeg.source.results.ctype = 'evoked';            % type of data: evoked=ERP, alternative e.g. 'trials'=single trial data
        matlabbatch{1}.spm.meeg.source.results.space = 1;
        matlabbatch{1}.spm.meeg.source.results.format = 'image';
        matlabbatch{1}.spm.meeg.source.results.smoothing = 8;               % FWHM smoothing kernel
        
        % run the job: nifti creation
        spm_jobman('run', matlabbatch)
    end

end

%% One sample Ttest on nifti images (only makes sense with multiple subjects!)
% -------------------------------------------------------------------------

suffix = '_f_1';                                                            % work with weird SPM suffixes
res_dir = fullfile(base_dir,'data', 'source_recon',out_fol);                % results directory

% subjects = [3,4];
% loop through time windows of interest
for w = 1:size(wois,1)
    
    subs_data = cell(1,numel(subjects));
    
    % loop through subjects and collect the path to nifti file 
    for s = 1:numel(subjects)    
        sub = sprintf('sub-%02d',subjects(s));
        subs_data{s} = fullfile(base_dir,'data', 'source_recon',sub,out_fol,        ...
                               [prefix sub '_SBL_'                      ...
                                num2str(w) '_t' num2str(wois(w,1)) '_' num2str(wois(w,2)) suffix '.nii']);
    end
    subs_data = cellstr(char(subs_data));
    
    save_dir = fullfile(res_dir,[prefix 't' num2str(wois(w,1)) '_' num2str(wois(w,2)) suffix]);
    if ~exist(save_dir,'dir')
        mkdir(save_dir)
    end
    
    clear matlabbatch   
    matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(save_dir);      % results directory  
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = subs_data;     % cellstr with nifti file paths
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;               % implicit masking
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};            % explicit masking: nifti file with spatial mask
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % run the job: specify statistical model on niftis
    spm_jobman('run', matlabbatch)
    
    % estimate model
    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr(fullfile(save_dir,'SPM.mat'));
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('serial',matlabbatch);    
end
