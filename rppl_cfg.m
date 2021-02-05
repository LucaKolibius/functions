% sets up the configuration before running the ripple detection
function tfg = rppl_cfg (data, ch)
tfg             = [];

% all channels at once (takes twice as long)
% tfg.thresCh     = data.label; % {'midHippR4'} or {'MiB1'};
% tfg.detectCh    = data.label(ch); %{'midHippR4'} or {'MiB1'};

% per channel
tfg.thresCh     = data.label(1);
tfg.detectCh    = data.label(1);

% tfg.param.bpfreq    = [80 120]; %lower and upper limit for bandpass filter % OLD
 tfg.param.bpfreq    = [80 140];
tfg.param.stageoi   = 1;
tfg.param.artfctPad = [-0.25 0.25]; %[prePad, postPad] (in s) additional padding around each
                                    % artifact, negative values = pre-artifact padding (default = [0 0])
tfg.param.thresType = 'channelwise';% Use thresholds corresponding to each channel, a mean threshold
                                    % across all channels or given fixed values. If channelwise is chosen
                                    % detectCh must be a subset of thresCh.
tfg.param.filtType  = 'fir'; % 'fir' or 'but' (default) - determines whether data processingis based on a FIR or butterworth filter
tfg.param.envType   = 'rms'; % 'env', 'rms' Root-mean-square level? (default), perform thresholding on rms signal or envelope
                             %(via spline interpolation of the peaks in the rectified signal) of bandpassed signal
tfg.param.envWin    = 0.02; % Window (in s) to use for envelope calculation (default 0.2 s)

tfg.criterion.len       = [0.038 0.5]; % [minLen, maxLen]: Minimal and maximal allowed duration of slow oscillation events
tfg.criterion.center    = 'mean'; % 'mean' (default) or 'median'
tfg.criterion.var       = 'centerSD';
%'centerSD' (default), 'centerMAD' (median absolute derivation), 'scaleCenter', 'scaleEnvSD' or 'percentile' or 'scaleFiltSD'
tfg.criterion.scale     = 2.5; % scalar representing the scaling factor applied to the 'criterionVar' parameter

tfg.paramOpt.smoothEnv   = 0.02; % If > 0 (default = 0), window (in s) to use for smoothing of envelope signal
tfg.paramOpt.upperCutoff = 9; % Scaling factor for an additional threshold, which discard any exceeding events to avoid artifacts (default = inf)
                              % Corresponding threshold is saved in cfg.thres.upperCutoff

tfg.paramOpt.minCycles       = 1; %Option to activate minimum number of required cycles per event
                                  % 0: Deactivated (default)
                                  % 1: Based on raw signal
                                  % 2: Based on bandpass filtered signal
tfg.paramOpt.minCyclesNum    = 3; % Scalar representing the minimum number of required cycles


%false positive rejection to eliminate IEDs
tfg.doFalseposRjct          = 1; % set to 1 if detected events are checked for their
                                 % frequency profile, i.e. a spectral peak with a
                                 % specific prominence within a specified frequency range
% tfg.falseposRjct.freqlim    = [75 125]; % [freqMin freqMax], frequency range event
%                                         % should have a spectral maximum
tfg.falseposRjct.freqlim    = [75 145]; % [freqMin freqMax], frequency range event
                                        % should have a spectral maximum
tfg.falseposRjct.timePad    = [-0.25 0.25]; % [timeMin timeMax], time padded around
                                            % events to calculate TFRs
tfg.falseposRjct.tfrFreq    = 65 : 2 : 155; % Frequencies for TFR calculation. Must include freqlim set above
tfg.falseposRjct.tfrTime    = -0.1 : 0.002 : 0.1; % Time range of TFR
tfg.falseposRjct.tfrWin     = ceil(0.1 * tfg.falseposRjct.tfrFreq) ./ tfg.falseposRjct.tfrFreq; % Length of time window for TFR calculation
tfg.falseposRjct.avgWin     = [-0.05 0.05]; % Time window used for averaging TFR, usual narrowly set around the event, i.e. t = 0 s
tfg.falseposRjct.prominence = 0.2; % Threshold the spectral peak amplitude has to exceed
