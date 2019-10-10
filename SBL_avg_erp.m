%% SBL: Average & Grand Average
% -------------------------------------------------------------------------
clc; clear; close all
addpath('C:\Users\Miro\Documents\Code\MATLAB\toolboxes\spm12')
spm('defaults', 'EEG');
project_dir = fullfile('D:','EEG','SCAN_class');

%% Single Subject Condition-Averages
% -------------------------------------------------------------------------

% prep data IDs
prefix = 'bletdhmar_';
SJs = [32,33];

% loop through subjects
for s = SJs
    
    subject = sprintf('sub-%02d',s);
    fprintf('Averaging %s \n',subject) 
    
    % prepare locations & data
    preproc_dir = fullfile(project_dir,'data',subject,'preprocessed');
    erp_dir = fullfile(project_dir,'analysis','erp',subject,'avg');
    if ~exist(erp_dir,'dir')
        mkdir(erp_dir)
    end
    D = spm_eeg_load(fullfile(preproc_dir, [prefix subject '_SBL.mat']));
       
    % copy data to separate location
    S = [];
    S.D = D;
    S.outfile = fullfile(erp_dir,D.fname);
    D = spm_eeg_copy(S);
               
    % get bad trials & catch trials & define as to delete
    load(fullfile(preproc_dir,['badtrials_' subject '_SBL.mat']));          % load badtrials (previously defined)
    events = D.events;                                                      % get event triggers from the meeg file
    events = cell2mat(cellfun(@(s)s.value,events,'uni',0));                 % save in array form
    catchtrials = find(events==33);                                         % get indices of catch trials with trigger number 33
    to_delete = sort(unique([badtrials,catchtrials]));                      % save indices in to_delete array
            
    % create average SEP (of all stimulation trials)
    % ---------------------------------------------------------------------
    %D = spm_eeg_load(fullfile(erp_dir,[prefix subject '_SBL.mat']));  % load initial file again 
    
    % define each trial as same condition for SEPs
    [new_conditions{1:size(D,3)}] = deal('SEP');
    D = conditions(D, ':', new_conditions);                                 % write conditions in meeg file with spm command   
    clear new_conditions                                                    % clear variable to avoid mistakes 
    
    % remove unwanted trials (careful, sorts trial order new!)
    D = D.badtrials(to_delete, 1);
    S           = [];
    S.D         = D;
    S.prefix    = 'r';
    D           = spm_eeg_remove_bad_trials(S);
    
    % average per conditions
    S           = [];
    S.robust    = 0;
    S.prefix    = 'SEP_';
    S.D         = D;
    D           = spm_eeg_average(S);
    
    % create average of HIGH & LOW 
    % ---------------------------------------------------------------------
    D = spm_eeg_load(fullfile(erp_dir,[prefix subject '_SBL.mat']));        % load initial file again 
    
    % define each trial as high/low according to trigger information (stored in events)   
    high_trials = find(events==12 | events== 22);
    low_trials = find(events==11 | events== 21);    
    [new_conditions{high_trials}] = deal('HIGH');
    [new_conditions{low_trials}] = deal('LOW');
    D = conditions(D, ':', new_conditions);                                 
    clear new_conditions
        
    % remove unwanted trials (careful, sorts trial order new!)
    D = D.badtrials(to_delete, 1);
    S           = [];
    S.D         = D;
    S.prefix    = 'r';
    D           = spm_eeg_remove_bad_trials(S);
    
    % average
    S           = [];
    S.robust    = 0;
    S.prefix    = 'HL_';
    S.D         = D;
    D           = spm_eeg_average(S);
  
    % create average of STANDARDS & DEVIANTS 
    % ---------------------------------------------------------------------
    D = spm_eeg_load(fullfile(erp_dir,[prefix subject '_SBL.mat']));        % load initial file again
    
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

    standard_trials = [stand_r1_hi_ind{:,s}, stand_r2_hi_ind{:,s}, ...
                       stand_r1_lo_ind{:,s}, stand_r2_lo_ind{:,s}];                   
    deviant_trials = [dev_r1_hi_ind{:,s}, dev_r2_hi_ind{:,s}, ...
                      dev_r1_lo_ind{:,s}, dev_r2_lo_ind{:,s}];   
                   
    [new_conditions{standard_trials}] = deal('STAND');
    [new_conditions{deviant_trials}] = deal('DEV');
    if size(new_conditions,1) < size(D,3)
        [new_conditions{size(D,3)}] = deal('Undefined');
    end
    D = conditions(D, ':', new_conditions);
    clear new_conditions
    
    % remove unwanted trials (careful, sorts trial order new!)
    undef = find(strcmp(D.conditions,'Undefined'));                         % additionally remove undefined trials (neither standard nor deviant)
    to_delete = sort(unique([to_delete, undef]));
    D = D.badtrials(to_delete, 1);
    S           = [];
    S.D         = D;
    S.prefix    = 'r';
    D           = spm_eeg_remove_bad_trials(S);
    
    % average
    S           = [];
    S.robust    = 0;
    S.prefix    = 'SD_';
    S.D         = D;
    D           = spm_eeg_average(S);    
       
end

%% Grand Averages
% -------------------------------------------------------------------------

% what avg file to do the GA on
GA_name = 'SD';
prefix = [GA_name '_rbletdhmar_'];
SJs = [6, 32, 33, 42];

% get file names for all subjects
data_files = [];
template_fid = {'nas', 'lpa', 'rpa'};
for s = 1:numel(SJs)
    subject = sprintf('sub-%02d',SJs(s));
    
    % Change fiducials labels to standard if unequal
    D = spm_eeg_load(fullfile(project_dir,'analysis','erp',subject,'avg',[prefix subject '_SBL.mat']));
    fid = D.fiducials;
    if ~isequal(fid.fid.label, template_fid)
      disp(subject); disp(fid.fid.label)
      fid.fid.label = template_fid;
      D = fiducials(D, fid);
      D.save;
    end
    
    data_files{s} = fullfile(project_dir,'analysis','erp',subject,'avg',[prefix subject '_SBL.mat']);    
end
data_files = char(data_files);

% GA
S = [];
S.D = data_files;
S.outfile = fullfile(project_dir,'analysis','erp',GA_name);

D = spm_eeg_grandmean(S);
