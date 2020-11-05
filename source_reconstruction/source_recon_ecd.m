%% Source reconstruction with SPM: Equivalent current dipole fitting (ECD)
% -------------------------------------------------------------------------
clc;clear;close

% Set Paths
addpath('C:\Users\Miro\Documents\Code\MATLAB\toolboxes\spm12')
base_dir = fullfile('E:','SCAN_class','neurocognitive_methods_class');

%% Grand Average ERP file
% -------------------------------------------------------------------------

fname = 'SEP.mat';                                                          % use SEP, SD
dname = 'P50';                                                              % use N20, P50, N140

% load GA
GA = fullfile(base_dir,'data','erp',fname);
D = spm_eeg_load(GA);

% make copy for fitting ECD
out_dir = fullfile(base_dir,'data','source_recon','ECD',dname);
if ~exist(out_dir,'dir')
    mkdir(out_dir) 
end 
S = [];
S.D = D;
S.outfile = fullfile(out_dir,fname);
D = spm_eeg_copy(S);

%% Forward Model
% -------------------------------------------------------------------------

fname = 'SD.mat';
dname = 'N140';

% load GA
GA = fullfile(base_dir,'data','source_recon','ECD',dname,fname);
D = spm_eeg_load(GA);

% remove any previous forward models
if isfield(D, 'inv')
    D = rmfield(D,'inv');
    D.save
end

% create forward model
clear matlabbatch
matlabbatch{1}.spm.meeg.source.headmodel.D = cellstr(GA);
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = 'SEP';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 1;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
spm_jobman('run', matlabbatch)

%% Optional: Inspect ERP of interest
% -------------------------------------------------------------------------

fname = 'SEP.mat';
dname = 'N20';

% load GA
GA = fullfile(base_dir,'data','source_recon','ECD',dname,fname);
D = spm_eeg_load(GA);

% convert to fieldtrip (convenient GUI)
data = fttimelock(D);

cfg               = [];
cfg.channel       = 'eeg';
close;
ft_topoplotER(cfg, data);

%% VB ECD (GUI)
% -------------------------------------------------------------------------

% RS2: 62 -12 14
% LS2: -62 -12 14

fname = 'SEP.mat';
dname = 'P50';

% load GA
GA = fullfile(base_dir,'data','source_recon','ECD',dname,fname);
D = spm_eeg_load(GA);

% % remove scaling
% D = rm_scale(D,1);
% D.save;

% start GUI: input = data and index of which forward model to base the fit on
spm_eeg_inv_vbecd_gui(D,1)

% display dipole in T1
D = spm_eeg_load(GA);
spm_eeg_inv_vbecd_disp('init',D,[],1)

%% Extract dipole waveforms
% -------------------------------------------------------------------------

fname = 'SD.mat';
dname = 'N140';

% load GA
GA = fullfile(base_dir,'data','erp',fname);
D = spm_eeg_load(GA);

% specify the locations and moments in dip.pnt and dip.ori respectively
dip.pnt = [41, -12, 53;
           36, -17, 44;
           46, -12, 16;
           -46, -12, 16];
dip.ori = [2.92, 1.67, 3.89 ;
           17.31, -14.25, -6.49;
           -28.31, 5.86, 4.24;
           12.69, 5.86, 4.24];
       
dip.label = {'s1_N20','s1_P50','rs2','ls2'};

S = [];
S.D = D;
S.dipoles = dip;

sD = spm_eeg_dipole_waveforms(S);

% save projections
out_dir = fullfile(base_dir,'data','source_recon','projections',dname);
if ~exist(out_dir,'dir')
    mkdir(out_dir) 
end 
save(fullfile(out_dir,'proj_dip_erp.mat'),'sD');