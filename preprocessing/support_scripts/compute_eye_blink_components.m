function compute_eye_blink_components(D, n_comps, new_fname)
% function that computes the priciple components in sensor space that
% correspond to eye-blinks
% D       - MEEG object which must hold information on eye-blinks (i.e. result of function detect_eye_blinks)
% n_comps - nr of principal components that should be calculated

if nargin < 3
    new_fname = 'ebf_conf';
end

% re-referance to common average
% S.D=D; 
% S.refchan='average';
% D_ebf1 = spm_eeg_reref_eeg(S);


% create epochs around all detected eye-blink events
S                         = [];
S.D                       = D;
S.timewin                 = [-200 400];
S.bc                      = 0;
S.trialdef.conditionlabel = 'blink';
S.trialdef.eventtype      = 'artefact';
S.trialdef.eventvalue     = {'eyeblink'};
% S.trialdef.conditionlabel = '7';
% S.trialdef.eventtype      = 'STATUS';
% S.trialdef.eventvalue     = 7;
S.reviewtrials            = 0;
S.save                    = 0;
S.epochinfo.padding       = 0;
D_ebf1                    = spm_eeg_epochs(S);


% mark bad trials as they might screw up the PCA
S                            = [];
S.D                          = D_ebf1;
S.badchanthresh              = 0.2;
S.methods.channels           = 'EEG';
S.methods.fun                = 'threshchan';
S.methods.settings.threshold = 800; % mV
D_ebf2                       = spm_eeg_artefact(S);


% compute the average of all eye-blink events
S        = [];
S.D      = D_ebf2;
S.robust = 0;
S.review = 0;
D_ebf3   = spm_eeg_average(S);  


% compute SVD to find first two principal components of avg. eye-blink in
% channel space
S         = [];
S.D       = D_ebf3;
S.method  = 'SVD';
S.timewin = [-200 400]; % in milliseconds!!! (for SPM12)!!!
S.ncomp   = n_comps;
D_ebf4    = spm_eeg_spatial_confounds_jh(S);

% copy data to simple file name
S         = [];
S.D       = D_ebf4;
S.outfile = new_fname;
spm_eeg_copy(S);

D_ebf1.delete();
D_ebf2.delete();
D_ebf3.delete();
