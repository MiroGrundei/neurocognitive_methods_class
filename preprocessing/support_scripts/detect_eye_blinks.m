function D_ebf = detect_eye_blinks(D, thresh, bipolar_chan_name, exg_chans)
% function to adjust threshold value for eye-blink correction. Works
% according to trial-and-error paired with visual inspection of the outcome
%
% D                  - Data object that was loaded
% thresh             - threshold for eye-blink detection (0 - 20 or even 50)
% bipolar_chan_name  - name of external electrode that will be recoded to
%                      hold bipolar information about eyeblinks (default = {'EXG8'})
% exg_chans          - names of external electrodes that are used to create
%                      bipolar channel. Cell-array of strings. (default = {'EXG3', 'EXG4'})

if nargin < 4
    exg_chans = {'EXG3', 'EXG4'};
end

if nargin < 3
    bipolar_chan_name = 'EXG8';
end

S = [];
S.D = D;
S.outfile = ['ebf_' D.fname];
D_ebf = spm_eeg_copy(S);

%turns 'bipolar_chan_name' channel into bipolar blink channel
D_ebf(D_ebf.indchannel(bipolar_chan_name),:) = ...
    D_ebf(D_ebf.indchannel(exg_chans{1}),:)-D_ebf(D_ebf.indchannel(exg_chans{2}),:); 

save(D_ebf);

%detect eyeblinks in channel specified by 'bipolar_chan_name'
S = [];
S.D = D_ebf;
S.eogchan = {bipolar_chan_name};
S.stdthresh= thresh;
S.overwrite = 1;
D_ebf = spm_eeg_detect_eyeblinks(S);
