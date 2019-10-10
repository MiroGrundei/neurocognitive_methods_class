%% SBL: Flexible Factorial GLM
% -------------------------------------------------------------------------
clc; close; clear all;
addpath('C:\Users\Miro\Documents\Code\MATLAB\toolboxes\spm12')
spm('defaults', 'EEG');
spm_jobman('initcfg');
project_dir = fullfile('D:','EEG','SCAN_class');

%% Prep files
% -------------------------------------------------------------------------

SJs = [32];
prefix = 'bletdhmar_';
img_name = 's12x12x0_condition_1.nii'; % s12x12x0_

% load condition information from saved indices 
index_dir = fullfile(project_dir, 'analysis', 'indices');
idx_fnames = {'stand_r1_hi_ind', ...
                   'stand_r1_lo_ind', ...
                   'stand_r2_hi_ind', ...
                   'stand_r2_lo_ind', ...
                   'dev_r1_hi_ind', ...
                   'dev_r1_lo_ind', ...
                   'dev_r2_hi_ind', ...
                   'dev_r2_lo_ind'};
for i = 1:numel(idx_fnames)
    load(fullfile(index_dir,[idx_fnames{i} '.mat']));
end 

%% Prep images & design matrix for all subjects
% -------------------------------------------------------------------------
 
X_sjs = cell(numel(SJs),1);
imgs_sjs = cell(numel(SJs),1);

for s = SJs
    
    subject = sprintf('sub-%02d',s);

    % prepare images
    data_dir = fullfile(project_dir,'data',subject,'preprocessed');
    img_file = fullfile(data_dir,['image_' prefix subject '_SBL'],img_name);
    imgs = cellstr(spm_select('expand',img_file));                          % spm command to expand the 3D image to multiple 2D filepaths
    
    % conditions of interest
    standards_r1 = [stand_r1_hi_ind{:,s}, stand_r1_lo_ind{:,s}];  
    standards_r2 = [stand_r2_hi_ind{:,s}, stand_r2_lo_ind{:,s}];    
    deviants_r1 = [dev_r1_hi_ind{:,s}, dev_r1_lo_ind{:,s}]; 
    deviants_r2 = [dev_r2_hi_ind{:,s}, dev_r2_lo_ind{:,s}];
    
    % specify design matrix X: 4000 trials by 2 factors
    X = zeros(size(imgs,1),2);

    % Factor 1: Standards & Deviants
    X([standards_r1,standards_r2],1) = 1;                                  
    X([deviants_r1,deviants_r2],1) = 2;                                     
    
    % Factor 2: Regime 1 & 2
    X([standards_r1,deviants_r1],2) = 1;                                  
    X([standards_r2,deviants_r2],2) = 2;  

    % throw out unspecified trials (w/o assigned condition)
    unspec = find( X(:,1)==0 & X(:,2)==0 );
    X(unspec,:) = [];
    imgs(unspec,:) = [];
    
    X_sjs{s} = X;                                                           % save for all subjects
    imgs_sjs(1:numel(imgs),s) = imgs;
    
end

%% Specify 1. Level GLM: Flexible Factorial
% -------------------------------------------------------------------------

glm_type = 'Regime_Type_Inter';

for s = SJs
    
    clear matlabbatch

    subject = sprintf('sub-%02d',s);
    fprintf('GLM %s \n',subject) 
    
    res_dir = fullfile(project_dir,'analysis','glm',subject,glm_type); 
    if ~exist(res_dir,'dir')
        mkdir(res_dir) 
    end 
    
    X = X_sjs{:,s};
    imgs = imgs_sjs(:,s);
    
    matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(res_dir);
    
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject.scans = imgs;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject.conds = X;          % mapping images to factors
                                                                                              % Factor and its levels defined in single column    
    % Factor 1
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'SD-Type';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 1;       % 0 = iid ; 1 => errors correlated within subject
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;   % 1 = unequal ; 0 => variance same for factor-levels
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;      % grand mean scaling
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;  
    
    % Factor 2
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Regime';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;       
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;   
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;      
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;    
    
    % Effects 
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [1 2];
    
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;               % implicit masking: 1 = masks outside of EEG layout
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};            % explicit masking: use a mask for specific time window or space of interest
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    spm_jobman('run', matlabbatch)                                          % run
    clear matlabbatch

    % Estimate
    disp(['Estimating GLM'])
    matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr(fullfile(res_dir,'SPM.mat'));
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    
    spm_jobman('serial',matlabbatch);
    clear matlabbatch  
    
end

disp('Done!')

%% 1. Level contrast
% -------------------------------------------------------------------------

% Design structure: [ SR1 | DR1 | SR2 | DR2 ]
con_name = 'D < S'; con_vec = [1 -1 1 -1];                                  % main effect of SD Type
% con_name = 'R1 < R2'; con_vec = [-1 -1 1 1];                                % main effect of regime
% con_name = 'SDtype_X_regime'; con_vec = [-1 1 1 -1];                        % SDtype X regime interaction

glm_type = 'Regime_Type_Inter';

for s = SJs

    clear matlabbatch

    subject = sprintf('sub-%02d',s);
    fprintf('GLM %s \n',subject) 
    
    res_dir = fullfile(project_dir,'analysis','glm',subject,glm_type);
                
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(fullfile(res_dir,'SPM.mat'));
        
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = con_name;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = con_vec;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
    spm_jobman('serial',matlabbatch);
    
end

%% Specify 2. Level GLM: One Sample T-Test (against H0 = no effect)
% -------------------------------------------------------------------------

img_name = 'con_0001.nii'; 
glm_type = 'Regime_Type_Inter';
res_dir = fullfile(project_dir,'analysis','glm');

imgs = cell(numel(SJs),1);
for s = 1:numel(SJs)   
    subject = sprintf('sub-%02d',SJs(s));    
    imgs{s} = fullfile(res_dir,subject,glm_type,img_name);     
end

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(res_dir);
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = imgs;

matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('serial',matlabbatch);           

% Estimate
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr(fullfile(res_dir,'SPM.mat'));
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('serial',matlabbatch);

