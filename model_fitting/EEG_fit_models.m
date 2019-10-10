%% Set-up Script for free energy approximation using VI-algorithm (spm_vi_glm.m) written by Dirk Ostwald
% directory prep
% -------------------------------------------------------------------------
clc;clear all;close all
addpath('
pdir        = ('/media/samg/DATA/Sam/SBL/data')                             ; % project directory
verbose     = 1                                                             ;

%  participant cohort parameters
% -------------------------------------------------------------------------
SJs         = [42]                                                          ; % participant indices
n_s         = length(SJs)                                                   ; % number of participants
n_t         = 359                                                           ; % number of peristimulus time bins

% data analysis parameters
% -------------------------------------------------------------------------
z_score     = 1                                                             ; % regressor zscoring flag
n_i         = 4                                                             ; % maximal number of VI algorithm iterations
delta       = 1e-4                                                          ; % variational free energy convergence criterion
alpha       = 1e-2                                                          ; % precision parameter of p(\beta) = N(\beta; 0, \alpha^{-1}I_p)
beta_lambda = 2e1                                                           ; % shape  parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
gama_lambda = 1e-1                                                          ; % scalar parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
dlab        = {'Cz', 'C6', 'Iz'}                                                       ; % electrode labels
n_d         = numel(dlab)                                                   ; % number of electrodes of interest

% data analysis array initializations
% -------------------------------------------------------------------------
plab        = {'TP'}                                                        ;
n_p         = numel(plab)                                                   ;
mlab        = {'PS', 'BS', 'CS'}                                            ; % model labels
mlab_full   = {'predictive_surprise','bayesian_surprise',                   ...
               'confidence_corrected_surprise'}                             ;
n_m         = 3                                                             ; % number of models of interest
Y_avg       = NaN(n_d,n_t,n_s)                                              ; % participant-specific trial-averaged dipole time-series
X           = cell(n_s,n_m)                                                 ; % participant-specific design matrices
n_trl       = NaN(n_s,1)                                                    ; % participant-specific number of trials upon bad and catch trial removal
lme         = NaN(n_t, n_m, n_p, n_d, n_s)                                  ; % log model evidence array

% participant iterations
% -------------------------------------------------------------------------
for s = 1:n_s

    % data preprocessing
    % ---------------------------------------------------------------------
    % participant-specific data filenames
    sdir  = fullfile(pdir, sprintf('sub-%02d', SJs(s)))                     ; % participant directory
    dfile = fullfile(sdir, 'preproc_anew',          ...
                     sprintf('bletdhmar_sub-%02d_SBL', SJs(s)))             ; % participant MEEG file
    bfile = fullfile(sdir, 'preproc_anew',                                  ...
                     sprintf('badtrials_sub-%02d_SBL.mat',SJs(s)))          ; % participant bad trials file

    
    % electrode data loading
    D           = spm_eeg_load(dfile);
    evts        = D.events;
    temp_events = ones(1, size(D.events,2));
    for i = 1:size(evts,2)
        temp_events(i) = evts{i}.value;
    end
    catch_trials = find(temp_events == 33);
    load(bfile, 'badtrials')                                                ; % bad trial indices
    time        = D.time                                                    ; % peristimulus time (identical over participants)
    Y           = D(:,:,:)                                                  ; % n_d x n_t x number of trials data array
    rmv         = unique(sort([catch_trials, badtrials]))                   ; % indices of trials to be removed
    n_trl(s)    = size(D,3) - length(rmv)                                   ; % valid trial number 
    Y(:,:,rmv)  = []                                                        ; % electrode data trial removal   
    

    
    % electrode iterations
    % ---------------------------------------------------------------------
    for d = 1:n_d
        d_idx = find(strcmp(D.chanlabels, dlab{d}))                         ; % find index of electrode
        
        % probability type iterations
        % -----------------------------------------------------------------
        for p = 1:n_p
            
            % model preprocessing
            % -------------------------------------------------------------   
            x = zeros(size(temp_events,2),5);            
            for m = 1:n_m                        
                mfile = fullfile(pdir, 'SCAN',  ...
                        sprintf('sub-%02d_tau_0.01_%s_CD.mat', ...
                                    SJs(s),plab{p}))                         ; % participant regressor file

                reg          = load(mfile)                                   ; % model regressor structure                                                                
                x(:,m)       = reg.(mlab_full{m})'                           ; % regressors of interest              
            end
            
            % regressor z scoring
            if z_score
                x = zscore(x);
            end
            
            % remove bad/catch-trials
            x(rmv,:) = [];                                              
            
            % model formulation & save DM
            n           = size(x,1)                                      ; % number of data points 
            DM          = cell(1,4);
            DM{1}       = ones(n,1)                                      ; % null model design matrix
            DM{2}       = [ones(n,1) x(:,1)]                             ; % predictive suprise model design matrix
            DM{3}       = [ones(n,1) x(:,2)]                             ; % bayesian suprise model design matrix
            DM{4}       = [ones(n,1) x(:,3)]                             ; % confidence corrected suprise model design matrix       
            X{s,p}      = DM;             

            % user information
            if verbose
                fprintf('Processing participant %d of %d, dipole %d of %d', s, n_s, d, n_d) 
                fprintf('\n')
            end
            
            % peri-stimulus time bin iterations
            % -------------------------------------------------------------
            for t = 1:n_t 
                
                % analysis model iterations 
                % ---------------------------------------------------------
                
                num_mdl = size(X{s,p},2);
                for mdl = 1:num_mdl
                
                    % spm_vb_glm analysis
                    glm             = []                                        ; % structure initialization
                    glm.X           = X{s,p}{mdl}                               ; % design matrix
                    glm.y           = squeeze(Y(d_idx,t,:))                     ; % data (constant over models)
                    glm.n           = size(glm.X,1)                             ; % number of data points
                    glm.p           = size(glm.X,2)                             ; % number of analysis model regression parameters
                    glm.mu_beta     = zeros(glm.p,1)                            ; % expectation parameter of p(\beta) = N(\beta;0_p,\alpha^{-1}I_p)
                    glm.alpha       = alpha                                     ; % precision parameter of p(\beta) = N(\beta; 0, \alpha^{-1}I_p)
                    glm.beta_lambda = beta_lambda                               ; % shape  parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
                    glm.gama_lambda = gama_lambda                               ; % scalar parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
                    glm.p           = size(glm.X,2)                             ; % number of beta parameters
                    glm.n_i         = n_i                                       ; % maximum number of iterations
                    glm.delta       = delta                                     ; % variational free energy convergence criterion
                    glm             = spm_vi_glm(glm)                           ; % estimation

                    % record converged variational free energy only
                    lme(t,mdl,p,d,s)    = glm.F_max; 
                    
                end
            end
        end
    end
end
disp('Done.')

%% Brief Plotting Example

% compute difference between models of interest and null-model
lme = squeeze(lme);
res = nan(n_t, n_m, n_d);
for i = 1:n_m % loop over models
    res(:,i, :) = lme(:,i+1, :) - lme(:,1, :);
end

clrmp = 'jet(256)'; % 'hot'; % 'jet(256)'; % jet(256)
figure(1);
n_mods = 3;
for d = 1:n_d
    % compute maximum relative F for axes
    mx = max(max(max(res)));
    if mx < 0
        mx = 1;
    end
    
    figure(1);            
    subplot(size(dlab,2),1,d)
    colormap(eval(clrmp))
    imagesc([res(:,1,d)'; res(:,2,d)'; res(:,3,d)'], [0 mx]); % mx(d)

    colorbar;
    set(gca,'YTick',[1:n_mods]);
    set(gca,'XTick',0:n_t/7:n_t); 
    set(gca,'YTickLabel',{'PS', 'BS', 'CS'}); 
    set(gca,'XTickLabel',[-100:100:600]);
    ylabel(dlab{d});
    
    sgtitle('Relative F (to null-model) across electrodes and peri-stimulus time - Sub-42')
    set(gca, 'FontSize', 11);  

end
%% Save
save(fullfile(pdir,'source_modeling','viglm','lme','lme_t-m-p-d-s.mat'),'lme')

