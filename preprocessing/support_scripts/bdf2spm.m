function D = bdf2spm(bdf_file,outfile,sj_loc_dir,ID)
% Convert biosemi .bdf file to data object for further processing with spm
% Includes sensor preparation for coregistration

% Read .bdf file
S = [];
S.dataset = bdf_file;
S.outfile = outfile;
D = spm_eeg_convert(S);

% Prepare sensor locations
sensor_file = fullfile(sj_loc_dir, [ID '_pos.mat']);
fid_file = fullfile(sj_loc_dir, [ID '_fid.mat']);
if exist(sensor_file,'file') && exist(fid_file,'file')
    S = [];
    S.D = D;
    S.sensfile = sensor_file;
    S.source = 'mat';
    S.headshapefile = fid_file;
    S.fidlabel = 'lpa nas rpa';
    S.task = 'loadeegsens';
    S.save = 1;
    D = spm_eeg_prep(S);
else % use default locations
    S = [];
    S.D = D;
    S.task = 'defaulteegsens';
    S.save = 1;
    D = spm_eeg_prep(S);

end

