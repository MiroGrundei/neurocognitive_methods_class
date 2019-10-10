function D = biosemi2spm(full_fname, outfile, default_locs, loc_folder)
% standard conversion of biosemi .bdf files to a data object for further
% processing with SPM (including the preparation for coregistration)
if nargin < 4
    loc_folder = '/mnt/storage/EEG_data/2012-12-10_EEG_Decision_Making_DMTS/loc_data';
end

if nargin < 3
    default_locs = [];
end

[fpath, fname, ext] = fileparts(full_fname);

subj_ID = char(regexp(fname, '\d+', 'match')); % a string
n = str2num(subj_ID);

S = [];
S.dataset = full_fname;
S.outfile = [outfile subj_ID];
D = spm_eeg_convert(S);

[junk file_mask junk] = fileparts(outfile);

loc_files = dir(fullfile(loc_folder, 'mat'));
if ismember([file_mask subj_ID 'pos.mat'], {loc_files.name}) && ...
   ismember([file_mask subj_ID 'fid.mat'], {loc_files.name}) && ...
   ~ismember(n,default_locs)
    S = [];
    S.D = D;
    S.sensfile = fullfile(loc_folder, 'mat', [file_mask subj_ID 'pos.mat']);
    S.source = 'mat';
    S.headshapefile = fullfile(loc_folder, 'mat', [file_mask subj_ID 'fid.mat']);
    S.fidlabel = 'lpa nas rpa';
    S.task = 'loadeegsens';
    S.save = 1;
    D = spm_eeg_prep(S);
% ... otherwise use default locations
else
    disp('>>>>>>>> DEFAULT POSITIONS <<<<<<<<<')
    S = [];
    S.D = D;
    S.task = 'defaulteegsens';
    S.save = 1;
    D = spm_eeg_prep(S);
end 