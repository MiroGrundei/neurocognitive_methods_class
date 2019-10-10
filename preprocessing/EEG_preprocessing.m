%% ========================================================================
%                   SCAN 2019 EEG PREPROCESSING
%% ========================================================================
% add paths containing spm12 and custom preprocessing scripts
addpath('~/programs/toolboxes/spm12')
addpath(genpath('/media/samg/DATA/Sam/SBL/data/SCAN/preprocessing'))
%addpath('/media/samg/DATA/Sam/SBL/code/spm_eeg')

% set data directory
data_dir        = fullfile('/media/samg/DATA/Sam/SBL/data/SCAN/preprocessing');

% project name
project         = 'SBL';

% load defaults
spm('defaults', 'EEG');

%% ========================================================================
%                       LOAD AND PREPARE RAW DATA
%% ========================================================================

% subjects to preprocess
SJs = [42];

for sj = 1:numel(SJs)    
    
    % prepare data
    % ---------------------------------------------------------------------
    subject     = sprintf('sub-%02d',SJs(sj));
    fprintf('Preprocessing %s \n',subject) % user info
        
    sj_dir      = fullfile(data_dir, subject, 'preprocessed');
    if ~exist(sj_dir,'dir')
        mkdir(sj_dir)
    end    
    loc_dir     = fullfile(data_dir,subject,'localization');
    bdf_file    = fullfile(data_dir,subject,'eeg',sprintf('%s_SBL.bdf',subject));
    out_file    = fullfile(sj_dir,sprintf('%s_SBL',subject));
    
    % convert location data (ZEBRIS) to .mat
    % ---------------------------------------------------------------------
    sfp2mat(loc_dir, sprintf('%s_SBL',subject))
    
    % convert data into SPM format & coregister electrode locations
    % ---------------------------------------------------------------------
    D           = bdf2spm(bdf_file,out_file,loc_dir,sprintf('%s_SBL',subject));
    
end

%% ========================================================================
%                               EDIT EVENTS
%% ========================================================================

SJs = [42];
prefix          = '';

for sj = 1:numel(SJs) 
    
    ID              = [subject '_' project];

    D               = spm_eeg_load(fullfile(sj_dir,[prefix ID '.mat']));

    % store current events
    evts            = D.events;
    new_evts        = D.events;
    vals            = [];

    % loop through all events
    for i           = 1:length(evts)
       if isempty(evts(i).value)
           evts(i).value = 999;            % change empty to some number, otherwise matlab is annoying
       end
       vals         = [vals,evts(i).value];% get values, save in list
       new_evts(i).value = [];             % extracted values get deleted
    end

    % find and classify different triggers
    end_trig        = find([evts.value] == 127);
    start_trig      = find([evts.value] == 126);
    time_trig       = find([evts.value] == 128);
    evt_trig        = find([evts.value] ==11|[evts.value]==12|[evts.value]==21| ...
                   [evts.value]==22|[evts.value]==33);

    % enter time trigger
    for i           = 1:length(time_trig)
       new_evts(time_trig(i)).value = vals(evt_trig(i));
    end

    % enter start trigger
    for i           = 1:length(start_trig)
       new_evts(start_trig(i)).value = 126;
    end

    % store recoded trigger in events and save
    D               = events(D, 1, new_evts);
    D.save;
end

%% ========================================================================
%                   MONTAGE > HPF > DOWNSAMPLING
%% ========================================================================

SJs               = [42];
for sj            = 1:numel(SJs)
    
    % load data
    % ---------------------------------------------------------------------
    subject       = sprintf('sub-%02d',SJs(sj));
    fprintf('Preprocessing %s \n',subject)                                 % user info

    sj_dir        = fullfile(data_dir, subject, 'preprocessed');
    
    D             = spm_eeg_load(fullfile(sj_dir, sprintf('%s_SBL.mat',subject)));
      
    % Re-reference to average montage 
    % ---------------------------------------------------------------------
    S             = [];
    S.D           = D;
    S.refchan     = 'average';
    S.prefix      = 'mar_'; 
    D             = spm_eeg_reref_eeg(S);

    % High-pass filter
    % ---------------------------------------------------------------------
    S             = [];
    S.D           = D;
    S.band        = 'high';
    S.freq        = 0.01;
    S.prefix      = 'h';
    D             = spm_eeg_filter(S);
       
    % Down-Sampling
    % ---------------------------------------------------------------------
    S             = [];
    S.D           = D;
    S.fsample_new = 512;
    S.prefix      = 'd';
    D             = spm_eeg_downsample(S);

end
disp('Done.')

%% ========================================================================
%                   EYE-BLINK DETECTION
%% ========================================================================

prefix          = '';
SJ              = 42;
subject         = sprintf('sub-%02d',SJ);

fprintf('Eye-blink detection %s \n',subject)   % user info

sj_dir          = fullfile(data_dir, subject, 'preprocessed');
D               = spm_eeg_load(fullfile(sj_dir, sprintf('%s%s_SBL.mat',prefix,subject)));
cd(sj_dir);

%% ========================================================================
%                   DETECT EYE-BLINKS 
%% ========================================================================
thresh          = 60;   % change HERE until appropriate value for subject is found
D_ebf           = detect_eye_blinks(D, thresh, 'EXG8', {'EXG3', 'EXG4'});
% execute this part twice: Once to compute first two components as a check 
% & finally only to compute the first component of eye-blink

%% ========================================================================
%        COMPUTE (first) COMPONENT OF EYE-BLINKS
%% ========================================================================
%The results of compute_eye_blink_components is automatically saved as 'ebf_config' in the subject folder.
num_components  = 1; 
compute_eye_blink_components(D_ebf, num_components);

%% eye-blink removal: loop through subjects
% -------------------------------------------------------------------------
threshs             = [60];
SJs                 = [42];

for sj              = 1:numel(SJs)

    prefix          = 'dhmar_';
    subject         = sprintf('sub-%02d',SJs(sj));
    fprintf('Eye-blink removal %s \n',subject)                                  % user info
    sj_dir          = fullfile(data_dir, subject, 'preprocessed');
    D               = spm_eeg_load(fullfile(sj_dir, sprintf('%s%s_SBL.mat',prefix,subject)));
    cd(sj_dir);
    
    % detect eye blinks
    % -------------------------------------------------------------------------
    thresh          = threshs(sj);
    D_ebf           = detect_eye_blinks(D, thresh, 'EXG8', {'EXG3', 'EXG4'});

    % principal components variance
    % -------------------------------------------------------------------------
    num_components  = 1; 
    compute_eye_blink_components(D_ebf, num_components); 
    
    close all
    
    % eye-blink removal
    % -------------------------------------------------------------------------
    try
        remove any spatial confounds file if present the meeg object
        S           = [];
        S.D         = D;
        S.method    = 'CLEAR';
        D           = spm_eeg_spatial_confounds(S);
    catch 
    end

    % add the spatial confound to the meeg object
    % -------------------------------------------------------------------------
    S               = [];
    S.D             = D;
    S.method        = 'SPMEEG';
    S.conffile      = 'ebf_conf.mat';
    D               = spm_eeg_spatial_confounds(S);

    % correct for the spatial confounds (Berg and Scherg)
    % -------------------------------------------------------------------------
    S               = [];
    S.D             = D;
    S.correction    = 'Berg';
    S.prefix        = 't';
    D               = spm_eeg_correct_sensor_data(S);
    
    close all
end

%% ========================================================================
%                EPOCHING, LOW-PASS FILTERING, BASELINE CORRECTION
%% ========================================================================

prefix              = 'tdhmar_';
SJs                 = [42];

for sj = 1:numel(SJs)
    
    subject         = sprintf('sub-%02d',SJs(sj));
    fprintf('Epoching %s \n',subject)       % user info
    sj_dir          = fullfile(data_dir, subject, 'preprocessed');
    D               = spm_eeg_load(fullfile(sj_dir, sprintf('%s%s_SBL.mat',prefix,subject)));
    
    % epoching
    % ---------------------------------------------------------------------
    S               = [];
    S.D             = D;
    S.bc            = 0;       % baseline correction: off
    S.timewin       = [-100 600];
    
    % specify events to be epoched
    evtlog          = [11,12,21,22,33];
    for j = 1:length(evtlog)
        S.trialdef(j).conditionlabel = '1';
        S.trialdef(j).eventtype      = 'STATUS';
        S.trialdef(j).eventvalue     = evtlog(j);
    end

    S.reviewtrials  = 0;
    S.prefix        = 'e';
    D               = spm_eeg_epochs(S);
    
    % LOW-PASS FILTERING
    % ---------------------------------------------------------------------
    S               = [];
    S.D             = D;
    S.band          = 'low';
    S.freq          = 45;
    S.prefix        = 'l';
    D               = spm_eeg_filter(S);
    
    % BASELINE CORRECTION
    % ---------------------------------------------------------------------
    S               = [];
    S.D             = D;
    S.timewin       = [-100 -5];
    S.prefix        = 'b';
    D               = spm_eeg_bc(S);
    
end

%% ========================================================================
%                   ARTEFACT DETECTION
%% ========================================================================
% Prepare and load file

prefix              = '';
SJ                  = 42;
subject             = sprintf('sub-%02d',SJ);

sj_dir              = fullfile(data_dir, subject, 'preprocessed');

D                   = spm_eeg_load(fullfile(sj_dir, sprintf('%s%s_SBL.mat',prefix,subject)));

%% VISUAL INSPECTION
% 'b': mark bad trials
% 's': save file
% 'escape': exit

thresh              = 100;
art_win             = [D.time(1) D.time(end)];
method              = 'threshchan';
dur                 = 8;          % duration in s for each plot
offset_val          = 50;   % y offsets between channels
D                   = eeg_inspect(D, thresh, art_win, method, dur, offset_val); % view epoched (for bad trial detection) or continuous data (for bad channel detection)  

%% save previously found bad trials
D                   = spm_eeg_load(fullfile(sj_dir, sprintf('%s%s_SBL.mat',prefix,subject)));
save(['badtrials_' subject '_' project '.mat'],'badtrials') % saves in current directory

%% Alternatively: Automatic Artefact Detection
% for our sub-42, it missed trials such as 102
S                            = [];
S.D                          = D;
S.badchanthresh              = 0.2;
S.methods.channels           = 'EEG';
%S.methods.fun                = 'jump';
S.methods.settings.threshold = 120; % muV
S.prefix                     = 'a';
S.methods.fun                = 'peak2peak';

D                            = spm_eeg_artefact(S);

%% ========================================================================
%                   OPTIONAL: BAD CHANNEL INTERPOLATION
%% ========================================================================

%addpath('~/programs/toolboxes/fieldtrip/template/neighbours/')
% Set bad channels using D.chanlabels in MEEG structure
chans = [35, 53]; 
D = D.badchannels(chans, 1); 

% perform linear interpolation
interpolate_bad_channels(D)

% remove bad channel labels for further processing
D = D.badchannels(chans, 0);
D.save();
