%% function to visually inspect single trial epoched data visually
%
% input
%   fileID  - full filename SPM EEG file
%   thresh  - threshold for artefact detection
%   art_win - window to check for artefacts
%   method  - method to check for artefacts ('threshchan' or 'jump' accordin to SPM)

function D = eeg_inspect(fileID, thresh, art_win, method, dur, offset_val)

%% Get data
if ~nargin || (~ischar(fileID) && ~isobject(fileID))
    [fileID, sts] = spm_select(1, 'mat', 'Select SPM M/EEG file');
    if ~sts; return; end
end
if ~isobject(fileID)
    D = spm_eeg_load(fileID);
else
    D = fileID;
end
if nargin < 2
    thresh = 80;
end
if nargin < 3
    art_win = [D.time(1) D.time(end)];
end
if nargin < 4
    method = 'threshchan';
end
if nargin < 5
    dur = 10;
end
if nargin < 6
    offset_val = 100;
end
allchan = D.indchantype('eeg');
nchan = numel(allchan);
badtrials = [];

ListenChar(2);

%% Restrict buttons for KbCheck
% KbName('UnifyKeyNames');
% Buttons = [ KbName('LeftArrow')     % Next view
%     KbName('RightArrow')    % Previous view
%     KbName('DownArrow')     % Next channel set
%     KbName('UpArrow')       % Previous channel set
%     KbName('Return')        % All channels
%     KbName('escape')        % Quit
%     KbName(',<')            % Backward skip n seconds
%     KbName('.>')            % Forward skip n seconds
%     KbName('b')             % Toggle bad trial
%     KbName('f')             % TF view
%     KbName('s')             % Save
%     KbName('t')             % Query trial number
%     KbName('e')             % Query epoch number
%     ];
% RestrictKeysForKbCheck(Buttons);


%% Open figure
main_fig = figure('name', 'EEG data inspector', 'units','normalized','outerposition',[0 0 1 1]);
% axes('position', [0 0.075 1 1])
axes('position', [0 0 1 1])
set(gca, 'visible', 'off')

try
    %% Epoched data
    if strcmp(D.type, 'single')
        
        i = 1;
        chansel = 1:nchan;
        t0 = find(D.time == 0);
        
        while 1
           % ListenChar(1);
            % get the data for given trial
            trial_data = squeeze(D(chansel,:,i));
            for c = 1:numel(chansel)
                trial_data(c,:) = trial_data(c,:) - trial_data(c,t0);
            end
            % find channels that have a signal above the threshold (abs value)
            if strcmpi(method, 'threshchan')
                above_thresh_ch = chansel((any(abs(trial_data(:,D.indsample(art_win(1)):D.indsample(art_win(2)))) > thresh, 2)));
            elseif strcmpi(method, 'jump')
                above_thresh_ch = chansel((any(abs(diff(trial_data(:,D.indsample(art_win(1)):D.indsample(art_win(2))), [], 2)) > thresh, 2)));
            else
                warning('Unknown artefact detection method: %s', method);
                above_thresh_ch = [];
            end
            
            % for plotting
            if length(chansel) == nchan
                %offset_val = 50; %sqrt(var(chunk_data(:)))*3;
            else
                offset_val = 100;
            end
            offsets = offset_val:offset_val:numel(chansel)*offset_val;
            offset_matrix = repmat(offsets, D.nsamples, 1);
            
            % plot EEG data
            figure(main_fig)
            plot(D.time*1000, trial_data+offset_matrix');
            hold on
            xlim([D.time(1)*1000 D.time(end)*1000])
            ylim([offsets(1)-offset_val*2 offsets(end)+offset_val*2])
            plot([0 0], ylim, 'k'); % trial onset
            hold off
            
            % mark good & bad channels, and unmarked channels that were above threshold
            badchan = intersect(chansel,D.badchannels);
            good_channels = setdiff(chansel, badchan);
            good_channels_below_thresh = setdiff(good_channels, above_thresh_ch);
            
            % Display channel labels
            % mark good channels below thresh in normal font
            if ~isempty(good_channels_below_thresh)
                idx = find(ismember(chansel,good_channels_below_thresh));
                text(repmat(D.time(10)*1000, 1, numel(offsets(idx))), offsets(idx), D.chanlabels(good_channels_below_thresh))
            end
            % mark good channels that are above threshold in this trial in bold, black font
            if ~isempty(above_thresh_ch)
                idx = find(ismember(chansel,above_thresh_ch));
                text(repmat(D.time(10)*1000, 1, numel(offsets(idx))), offsets(idx), D.chanlabels(above_thresh_ch), 'color', 'k', 'Fontweight', 'bold', 'FontSize', 12)
            end
            % mark bad channels in red and large font
            if ~isempty(badchan)
                idx = find(ismember(chansel,badchan));
                text(repmat(D.time(10)*1000, 1, numel(offsets(idx))), offsets(idx), D.chanlabels(badchan), 'color', 'r', 'Fontweight', 'bold', 'FontSize', 15)
            end
            
            % Display trial number
            text(D.time(round(end/2))*1000, offsets(round(end/2)), sprintf('trial #%d',i), 'FontSize', 25, 'color', [.5 .5 .5])
            
            % Display bad trial label
            if D.badtrials(i) %~isempty(intersect(above_thresh_ch, good_channels))
                text(D.time(round(end/2))*1000, offsets(round(end/2)-2), 'bad trial', 'FontSize', 25, 'color', [1 0 0])
            end
            
            drawnow
            
            
            %% Read keyboard
            
            [~, keyCode, ~] = KbStrokeWait();
            keyName = KbName(keyCode);
            
            switch keyName
                
                % Quit
                case 'ESCAPE'
                    save(D);
                    close(main_fig)
                    ListenChar(1);
                    return
                    
                    % Next trial
                case 'RightArrow'
                    i = i + 1;
                    if i > D.ntrials
                        i = 1;
                    end;
                    
                    % Previous trial
                case 'LeftArrow'
                    i = i - 1;
                    if i < 1
                        i = D.ntrials;
                    end
                    
                    % Next channel set
                case 'DownArrow'
                    n = 32; % Number of channels displayed at once (should be factor of 64!)
                    if isempty(setdiff(allchan,chansel))
                        chansel = (nchan-(n-1)):nchan;
                    else
                        chansel = (chansel(1)-n):(chansel(1)-1);
                        if any(chansel < 1)
                            chansel = (nchan-(n-1)):nchan;
                        end
                    end
                    
                    % Previous channel set
                case 'UpArrow'
                    n = 16; % Number of channels displayed at once (should be factor of 64!)
                    if isempty(setdiff(allchan,chansel))
                        chansel = 1:n;
                    else
                        chansel = (chansel(end)+1):(chansel(end)+n);
                        if any(chansel > nchan)
                            chansel = 1:n;
                        end
                    end
                    
                    % All channels
                case 'Return'
                    chansel = D.indchantype('eeg');
                    
                    % Toggle bad trial
                case 'b'
                    if D.badtrials(i)
                        D = D.badtrials(i,0);
                    else
                        D = D.badtrials(i,1);
                    end
                    
                    % query trial number
                case 't'
                    ListenChar(1);
                    i = input(sprintf('Trial number (1-%s): ',num2str(D.ntrials)));
                    ListenChar(2);
                    if i < 1
                        i = 1;
                    elseif i > D.ntrials
                        i = D.ntrials;
                    end
                    
                    % save the data
                case 's'
                    save(D);
                    disp('Saved dataset!')
                    
                    % compute TF of trial with standard specs: morlet, 25 ms resolution
                case 'f'
                    S = [];
                    S.subsample   = 25;
                    S.ncycles     = 7;
                    S.frequencies = [4:2:48];
                    
                    res = spm_eeg_specest_morlet(S, trial_data, D.time);
                    
                    pow = res.fourier.*conj(res.fourier);
                    
                    inds   = find(res.time <= -2.25, 1, 'last'):find(res.time >= -1.5, 1, 'first');
                    xbase  = mean(pow(:,:,inds),3);
                    bc_pow = 100*((pow./repmat(xbase,[1 1 size(pow,3)]) - 1));
                    
                    chsel = inputdlg('Which channel?', 'Channel selection');
                    
                    if isempty(chsel), continue, end
                    
                    figure('name', sprintf('trial #%d (%s) - cond: %s', i, chsel{1}, D.conditions{i}));
                    
                    imagesc(res.time, res.freq, squeeze(pow(D.indchannel(chsel), :, :))); %, [0 1000])
                    colorbar;
                    xlabel('time')
                    ylabel('freq')
                    axis xy
                    drawnow;
                    
                    continue;
            end
            
        end
        
        ListenChar(1);
        
    elseif strcmp(D.type, 'continuous')
        
        % loop over the data in chunks
%        dur = 10; % secs
        i = 0; % in seconds
        e = 1;
        chansel = 1:nchan;
        
        % get epochs (runs etc.)
        evts = D.events;
        for ev = 1:numel(evts)
            if isempty(evts(ev).value)
                evts(ev).value = 888;
            end
        end
        evt_t = [evts.time];
        evt_vals = [evts.value];
        eps = find(strcmp({evts.type},'Epoch'));
        eps = [eps numel(evts)];
        eps_t = [evts(eps).time];
        eps_t = [evts.time];
        nepochs = length(eps);
        
        while 1
            
            chunk_samples = D.indsample(i):D.indsample(i+dur);
                        
            % get the data for given trial
            chunk_data = squeeze(D(chansel,chunk_samples));
                        
            % find channels that have a signal above the threshold (abs value)
            if strcmpi(method, 'threshchan')
                above_thresh_ch = chansel(any(abs(chunk_data) > thresh, 2));
            elseif strcmpi(method, 'jump')
                above_thresh_ch = chansel(any(abs(diff(chunk_data, [], 2)) > thresh, 2));
            else
                warning('Unknown artefact detection method: %s', method);
                above_thresh_ch = [];
            end
                        
            % for plotting
            if length(chansel) == nchan
%                 offset_val = 100;
                %offset_val = sqrt(var(chunk_data(:)))*5;
            else
                offset_val = 100;
            end
            offsets = offset_val:offset_val:numel(chansel)*offset_val;
            
            offset_matrix = repmat(offsets, numel(chunk_samples), 1);
            
            % ========== plot EEG data ===============
                        
            % get the right time
            t = D.time(chunk_samples);

            % get the events that happen in that chunk
            evt_idx = find((evt_t >= t(1)) & (evt_t <= t(end)));
            evts_of_interest = find([evts.value]==11| ...
                        [evts.value]==12| ...
                        [evts.value]==21| ...
                        [evts.value]==22| ...
                        [evts.value]==33);
            eoi_idx_label = find(ismember(evts_of_interest, evt_idx));
            
            chunk_evt_vals = [];
            for c = 1:numel(evt_idx)
                chunk_evt_vals = [chunk_evt_vals evts(evt_idx(c)).value]; 
            end
            chunk_evt_idx = evt_idx(find(ismember(chunk_evt_vals,[11,12,21,22,33])));
            
            plot(t, chunk_data+offset_matrix', 'LineWidth', 0.5);
            hold on
            xlim([t(1) t(end)])
            ylim([offsets(1)-offset_val*2 offsets(end)+offset_val*2])
            plot([evt_t(chunk_evt_idx)-0.1; evt_t(chunk_evt_idx)-0.1], repmat(ylim, numel(chunk_evt_idx), 1)', 'k')
            hold off

            for j=1:numel(chunk_evt_idx)
                %text(evt_t(evt_idx(j)), 0, num2str(evt_vals(evt_idx(j))), 'Fontsize', 14)
                text(evt_t(chunk_evt_idx(j)), 0, num2str(eoi_idx_label(j)), 'Fontsize', 14)
            end
            
            % mark good & bad channels, and channels that were above threshold
            good_channels_below_thresh = setdiff(chansel, above_thresh_ch);
            
            % mark good channels below thresh in normal font
            if ~isempty(good_channels_below_thresh)
                idx = find(ismember(chansel,good_channels_below_thresh));
                text(repmat(t(10), 1, numel(offsets(idx))), offsets(idx), D.chanlabels(good_channels_below_thresh),'FontSize', 12)
            end
            
            % mark good channels that are above threshold in this trial in bold, red font
            if ~isempty(above_thresh_ch)
                idx = find(ismember(chansel,above_thresh_ch));
                text(repmat(t(10), 1, numel(offsets(idx))), offsets(idx), D.chanlabels(above_thresh_ch), 'color', 'r', 'Fontweight', 'bold', 'FontSize', 12)
            end
            
%             ListenChar(1)
            yl = ylim;
            text(t(500), yl(2)*.95, sprintf('Time: %s - %s sec', num2str(round(i)), num2str(round(i+dur))), 'FontSize', 24, 'color', [0 0 0])
            text(t(500), yl(2)*.9, sprintf('Epoch: %d', e), 'FontSize', 24, 'color', [0 0 0])
%             if ~isempty(above_thresh_ch)
%                 text(t(round(end/2)), offsets(round(end/2)-2), 'bad chunk', 'FontSize', 25, 'color', [1 0 0])
%             end
            
            drawnow
            
            %% Read keyboard
            
            [~, keyCode, ~] = KbStrokeWait();
            keyName = KbName(keyCode);
            
            switch keyName
                
                % Quit
                case 'ESCAPE'
                    % save(D);
                    close(main_fig)
                    ListenChar(1);
                    return
                    
                    % Next trial
                case 'RightArrow'
                    i = i + dur;
                    if i+dur > D.time(end)
                        i = 0;
                    end;
                    tmp = sort([eps_t i]);
                    idx = find(tmp == i);
                    if length(idx) > 1
                        idx = idx(end);
                    end
                    if idx == 1
                        e = 1;
                    else
                        ep_t = tmp(idx-1);
                        e = find(eps_t == ep_t);
                    end
                    
                    % Previous trial
                case 'LeftArrow'
                    i = i - dur;
                    if i < 0
                        i = D.time(end)-dur;
                    end
                    tmp = sort([eps_t i]);
                    idx = find(tmp == i);
                    if length(idx) > 1
                        idx = idx(end);
                    end
                    if idx == 1
                        e = 1;
                    else
                        ep_t = tmp(idx-1);
                        e = find(eps_t == ep_t);
                    end
                    
                % Forward skip
                case '.>'
                    i = i + dur*10;
                    if i+dur > D.time(end)
                        i = 0;
                    end;
                    tmp = sort([eps_t i]);
                    idx = find(tmp == i);
                    if length(idx) > 1
                        idx = idx(end);
                    end
                    if idx == 1
                        e = 1;
                    else
                        ep_t = tmp(idx-1);
                        e = find(eps_t == ep_t);
                    end
                    
                    % Backward skip
                case ',<'
                    i = i - dur*10;
                    if i < 0
                        i = D.time(end)-dur;
                    end
                    tmp = sort([eps_t i]);
                    idx = find(tmp == i);
                    if length(idx) > 1
                        idx = idx(end);
                    end
                    if idx == 1
                        e = 1;
                    else
                        ep_t = tmp(idx-1);
                        e = find(eps_t == ep_t);
                    end
                    
                    % Next channel set
                case 'DownArrow'
                    n = 32; % Number of channels displayed at once (should be factor of 64!)
                    if isempty(setdiff(allchan,chansel))
                        chansel = (nchan-(n-1)):nchan;
                    else
                        chansel = (chansel(1)-n):(chansel(1)-1);
                        if any(chansel < 1)
                            chansel = (nchan-(n-1)):nchan;
                        end
                    end
                    
                    % Previous channel set
                case 'UpArrow'
                    n = 16; % Number of channels displayed at once (should be factor of 64!)
                    if isempty(setdiff(allchan,chansel))
                        chansel = 1:n;
                    else
                        chansel = (chansel(end)+1):(chansel(end)+n);
                        if any(chansel > nchan)
                            chansel = 1:n;
                        end
                    end
                    
                    % All channels
                case 'Return'
                    chansel = D.indchantype('eeg');
                                        
                    % query epoch number
                case 'e'
                    ListenChar(1);
                    e = input(sprintf('Epoch (Block) number (1-%s): ',num2str(nepochs)));
                    ListenChar(2);
                    if e < 1
                        e = 1;
                    elseif e > nepochs
                        e = nepochs;
                    end
                    
                    i = evts(eps(e)).time-dur;
                    
                    % save badtrial epoch number
                case 'b'    
                    ListenChar(1);
                    evt_idx_inpt = input('Trigger Index: ');                  
                    ListenChar(2);
                    badtrials = [badtrials, evt_idx_inpt];
                    save(fullfile(fileparts(fileID),'badtrial_idx.mat'),'badtrials');
                    
            end
            
        end
        
        ListenChar(1);
        
    else
        
        % can't handle the data type
        
        fprintf('***** Sorry, function does not support datatype: %s', D.type);
        
        ListenChar(1);
        
    end
    
catch e
    ListenChar(1);
    disp(e)
end
