%% SBL: Plot ERP
% -------------------------------------------------------------------------
clc; clear; close all
addpath('C:\Users\Miro\Documents\Code\MATLAB\toolboxes\spm12')
spm('Defaults','EEG')
project_dir = fullfile('D:','EEG','SCAN_class');

%% Plot Subject Average ERP
% -------------------------------------------------------------------------

% prepare data & plot info
subject = 'sub-42';
condition = 'SD';
channel = 'CP6';                                                            % channel name or 'all'
prefix = 'rbletdhmar_';
erp_dir = fullfile(project_dir,'analysis','erp');

% load data
data = fullfile(erp_dir,subject,'avg',[condition '_' prefix subject '_SBL.mat']);
D = spm_eeg_load(data);
data_GA = fullfile(erp_dir,[condition '.mat']);                             % grand average data
D_GA = spm_eeg_load(data_GA);

% plot
if strcmp(channel,'all')                                                    % overview of all channels   
    for i = 1:64
        for cond = 1:size(D,3)
            subplot(8,8,i)
            plot(D.time,D(i,:,cond),'linewidth', 1);
            hold on
        end
        title(sprintf('%s',D.chanlabels{i}))
        hold off
    end
    legend(D.conditions)
else                                                                        % specific channel     
    nChan = find(strcmp(D.chanlabels,channel));
    figure(1); 
    for cond = 1:size(D,3)
        plot(D.time,D(nChan,:,cond),'linewidth', 2);
        hold on           
    end
    title(sprintf('Average %s %s at electrode %s',condition,subject,channel))
    legend(D.conditions)
    hold off    
    
    % plot grand average
    figure(2); 
    for cond = 1:size(D,3)
        plot(D_GA.time,D_GA(nChan,:,cond),'linewidth', 2);
        hold on           
    end
    title(sprintf('Grand average %s at electrode %s',condition,channel))
    legend(D.conditions)
    hold off
end
 