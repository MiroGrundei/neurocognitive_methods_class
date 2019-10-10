% function to interpolate bad channels according to biosemi neighborhood relationships of channels

function D = interpolate_bad_channels(D)

% prepare neighbourhood structure of electrodes
cfg = [];
cfg.method = 'template';
cfg.template = 'biosemi64_neighb.mat';

elecs = ft_prepare_neighbours(cfg);

% create a copy of the dataset
% S         = [];
% S.D       = D;
% S.outfile = ['intp' D.fname]; % SPM12
% S.newname = ['intp' D.fname]; % SPM8
% D_intp    = spm_eeg_copy(S);

for bad_idx = D.badchannels
    
    neighbours = elecs(strcmpi({elecs.label}, D.chanlabels(bad_idx))).neighblabel;
    
    D(bad_idx,:,:) = mean(D(D.indchannel(neighbours),:,:),1);
    
end