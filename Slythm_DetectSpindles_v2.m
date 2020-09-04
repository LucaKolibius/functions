%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slythm_DetectSOs
% by T.O. Bergmann, modified by H.-V.V Ngo & B. Staresina
%
% Detects discrete slow oscillation events in EEG data
% - Requires Fieldtrip format
% - Requires hypnogram within EEG data
% - If artifacts ('Channel x Trials' cell with each cell containing a [numArt x 2] 
%   array of artifact start & end) are included in inData, they will be
%   considered to discard unfit event candidates.
%
% Usage: outEvt = Slythm_DetectSpindles_v2(cfg, inData)
%
% Configuration parameters
% cfg.thresCh       = (N x 1) string listing the N channels to determine the thresholds for
% cfg.detectCh      = (M x 1) string listing the M (< N) channels to detect SO in (based on the previously determined thresholds)

% cfg.param.stageoi     = (N_SlSt x 1) vector containing sleep stages of interest for detection
% cfg.param.bpfreq      = [freqMin, freqMax]: Lower and upper limit for bandpass filter
% cfg.param.filtType    = 'fir' or 'but' (default) - determines whether data processingis based on a FIR or butterworth filter
% cfg.param.thresType   = 'channelwise' (default), 'average' or 'fixed' (not implemented)
%                         Use thresholds corresponding to each channel, a mean threshold 
%                         across all channels or given fixed values. If channelwise is chosen
%                         detectCh must be a subset of thresCh.
% cfg.param.envType     = 'env', 'rms' (default), perform thresholding on rms signal or envelope (via spline interpolation of the peaks in the rectified signal) of bandpassed signal
% cfg.param.envWin      = Window (in s) to use for envelope calculation (default 0.2 s)
% cfg.param.artfctPad   = [prePad, postPad] (in s) additional padding around each 
%                         artifact, negative values = pre-artifact padding (default = [0 0])
%
% cfg.criterion.len      = [minLen, maxLen]: Minimal and maximal allowed duration of slow oscillation events
% cfg.criterion.center   = 'mean' (default) or 'median'
% cfg.criterion.var      = 'centerSD' (default), 'centerMAD' (median absolute derivation), 'scaleCenter', 'scaleEnvSD' or 'percentile' or 'scaleFiltSD'
% cfg.criterion.scale    = scalar representing the scaling factor applied to the 'criterionVar' parameter
% cfg.criterion.padding  = [prePad, postPad] (in s) additional padding around
%                         each event candidate, negative value = pre-event 
%                         padding (default = [0 0])

% cfg.paramOpt.smoothEnv    = If > 0 (default = 0), window (in s) to use for smoothing of envelope signal
% cfg.paramOpt.upperCutoff  = Scaling factor for an additional threshold, which discard any exceeding events to avoid artifacts (default = inf)
%                             Corresponding threshold is saved in cfg.thres.upperCutoff
% cfg.paramOpt.mergeEvts    = if > 0, maximal gap (in s) until which close events are merged
% cfg.paramOpt.scndThres    = if > 0 (default = 0), introduce a second
%                             (higher) threshold criterion threshold saved
%                             in cfg.thres.second
% cfg.paramOpt.minCycles    = Option to activate minimum number of required cycles per event
%                             0: Deactivated (default)
%                             1: Based on raw signal
%                             2: Based on bandpass filtered signal
% cfg.paramOpt.minCyclesNum = Scalar representing the minimum number of required cycles
%
% cfg.doFalseposRjct  = set to 1 if detected events are checked for their
%                         frequency profile, i.e. a spectral peak with a
%                         specific prominence within a specified frequency range
% cfg.falseposRjct.freqlim    = [freqMin freqMax], frequency range event
%                                 should have a spectral maximum
% cfg.falseposRjct.timePad    = [timeMin timeMax], time padded around
%                                 events to calculate TFRs
% cfg.falseposRjct.tfrFreq    = Frequencies for TFR calculation. Must
%                                 include freqlim set above
% cfg.falseposRjct.tfrTime    = Time range of TFR
% cfg.falseposRjct.tfrWin     = Length of time window for TFR calculation
% cfg.falseposRjct.avgWin     = Time window used for averaging TFR, usual
%                                 narrowly set around the event, i.e. t = 0 s
% cfg.falseposRjct.prominence = Threshold the spectral peak amplitude has to exceed
%
% example input structure for two channel analyis:
% inData.label:       {'Fz'; 'Cz'}
% inData.time:        {[2×37320000 double]}
% inData.trial:       {[2×37320000 double]}
% inData.staging:     {[1×37320000 int8]}
% inData.fsample:     1000
% inData.artifacts:   {[46020×2 double]; [50020×2 double]}
% inData.sampleinfo:  [1 37320000]
%
%
% last update 18-11-09 by HVN
% To do:
% [ ] Add description of false-positive rejection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cfg = hvn_slythm_detectSpindle_v2(cfg, inData)

%% Check numbers of arguments
if nargin ~= 2
    error('Wrong number of arguments');
end


%% Check input channels
if ~isfield(cfg,'thresCh') || ~isfield(cfg,'detectCh')
    error('No thresCh or detectCh specified');
end

if ismember('all', cfg.thresCh)
    cfg.thresCh = inData.label;
end

if ismember('all', cfg.detectCh)
    cfg.detectCh = inData.label;
end

if strcmp(cfg.param.thresType, 'channelwise') %% Check if detectCh is a subset of thresCh
    if sum(ismember(cfg.thresCh, cfg.detectCh)) ~= size(cfg.detectCh)
        error('Channelwise detection requires detectCh to be a subset of thresCh');
    end
end


%% Housekeeping
%--- Check if scoring is available
if ~isfield(inData, 'staging')
    error('No scoring available');
end

%--- Check if cfg.criterion.center exists, otherwise set to 'mean'
if ~isfield(cfg.criterion,'center')
    cfg.criterion.center = 'mean';
end

%--- Check if filter type is specified
if isempty(cfg.param.filtType)
    cfg.param.filtType = 'but';
end

%--- Check of envelope parameters exists
if ~isfield(cfg.param,'envType')
    cfg.param.envType = 'rms';
end

if ~isfield(cfg.param,'envWin')
    cfg.param.envWin = 0.2;
end

%--- Check smoothing parameter
if ~isfield(cfg.paramOpt, 'smoothEnv')
    cfg.paramOpt.smoothEnv = 0;
end

%--- Check if a second threshold criterion is applied
if ~isfield(cfg.paramOpt, 'scndThres')
    cfg.paramOpt.scndThres = 0;
end

%--- Check padding
if ~isfield(cfg.criterion,'padding')
    cfg.criterion.padding = [0 0];
end

if ~isfield(cfg.param,'artfctPad')
    cfg.param.artfctPad = [0 0];
end

%--- Check if upperCutoff is existing
if ~isfield(cfg.paramOpt, 'upperCutoff')
    cfg.paramOpt.upperCutoff = inf;
end

%--- Check if event candidates shall be merged
if ~isfield(cfg.paramOpt,'mergeEvts')
    cfg.paramOpt.mergeEvts = 0;
end

%--- Check minCycle option
if ~isfield(cfg.paramOpt,'minCycles')
    cfg.paramOpt.minCycles      = 0;
    cfg.paramOpt.minCyclesNum   = 0;
elseif cfg.paramOpt.minCycles > 0 && ~isfield(cfg.paramOpt, 'minCyclesNum')
    error('Specification of minCyclesNum is missing');
end

if ~isfield(cfg,'doFalseposRjct'); cfg.doFalseposRjct = 0; end


%% Prepare output structure
cfg.thres       = [];
cfg.evtIndiv    = struct;
cfg.evtSummary  = struct;


%% Important variables
fsample         = round(inData.fsample);
numTrl          = size(inData.trial,2);
numThresCh      = size(cfg.thresCh,1);
numDetectCh     = size(cfg.detectCh,1);
artfctPad       = round(cfg.param.artfctPad * fsample);
criterionPad    = round(cfg.criterion.padding * fsample);


%% Prepare binary vector representing samples within desired sleep stages
bnryStage = cellfun(@(x) ismember(x, cfg.param.stageoi),inData.staging,'UniformOutput', 0);


%% Prepare artfctFilter
bnryArtfree = cellfun(@(x) ones(numThresCh,size(x,2),'logical'),inData.trial,'UniformOutput',0);

if isfield(inData, 'artifacts')
    for iTrl = 1 : numTrl
        for jCh = 1 : numThresCh
            currCh = ismember(inData.label,cfg.thresCh{jCh});
            
            % @ LDK
            if isempty(inData.artifacts{currCh,iTrl})
                continue
            end
            
            tmpArt = inData.artifacts{currCh,iTrl} + artfctPad;
            tmpArt(tmpArt(:,1) < 1,1)                         = 1;                              %% Ensure padding does not...
            tmpArt(tmpArt(:,2) > inData.sampleinfo(iTrl,2),2) = inData.sampleinfo(iTrl,2);      %% exceed data range
            
            for iArt = 1 : size(tmpArt,1)
                bnryArtfree{iTrl}(jCh,tmpArt(iArt,1):tmpArt(iArt,2)) = 0;
            end
            
            clear tmpArt
        end
    end   
end


%% Prepare data to determine thresholds
tmpThres = cell(1,numTrl);

for iTrl = 1 : numTrl                                                   % loop over trials
    for jCh = 1 : numThresCh                                            % loop over channels
        fprintf('Filter channel %d/%d (%s)\n', jCh, numThresCh, cfg.thresCh{jCh});
        
        currCh = ismember(inData.label,cfg.thresCh{jCh});
        switch cfg.param.filtType
            case 'fir'
                tmpThres{iTrl}(jCh,:) = ft_preproc_bandpassfilter(inData.trial{iTrl}(currCh,:), fsample, cfg.param.bpfreq, 3*fix(fsample/cfg.param.bpfreq(1))+1, 'fir', 'twopass');
            case 'but'
                tmpHP = ft_preproc_highpassfilter(inData.trial{iTrl}(currCh,:), fsample, cfg.param.bpfreq(1), 5, 'but', 'twopass','reduce');
                tmpThres{iTrl}(jCh,:) = ft_preproc_lowpassfilter(tmpHP, fsample, cfg.param.bpfreq(2), 5, 'but', 'twopass','reduce');
                
                clear tmpHP
        end
    end
end


%% calculate envelope of bandpass filtered signal
tmpThresEnv = [];
switch cfg.param.envType
    case 'rms'
        tmpThresEnv.trial = cellfun(@(x) envelope(x',round(cfg.param.envWin * fsample),'rms')',tmpThres,'UniformOutput',0);
    case 'env'
        tmpThresEnv.trial = cellfun(@(x) envelope(x',round(cfg.param.envWin * fsample),'analytic')', tmpThres,'UniformOutput',0);
end

if cfg.paramOpt.smoothEnv > 0 %% Smooth envelope if requested
    if ~verLessThan('matlab','9.2')
        tmpThresEnv.trial = cellfun(@(x) smoothdata(x,2,'movmean',round(cfg.paramOpt.smoothEnv * fsample)),tmpThresEnv.trial,'UniformOutput',0);
    else
        for iTrl = 1 : numTrl
            for jCh = 1 : numThresCh
                tmpThresEnv.trial{iTrl}(jCh,:) = smooth(tmpThresEnv.trial{1,iTrl}(jCh,:),round(cfg.paramOpt.smoothEnv * fsample));
            end
        end
    end
end


%% determine amplitude threshold
% gather data for variance threshold calculation
vardata = cell(numThresCh,1);
for iTrl = 1 : numTrl                           % loop over trials
    if strcmp(cfg.criterion.var,'scaleFiltSD')
        for jCh = 1 : numThresCh
            vardata{jCh} = [vardata{jCh}, tmpThres{iTrl}(jCh,all([bnryStage{iTrl}; bnryArtfree{iTrl}(jCh,:)]))];
        end
    else
        for jCh = 1: numThresCh
            vardata{jCh} = [vardata{jCh}, tmpThresEnv.trial{iTrl}(jCh,all([bnryStage{iTrl}; bnryArtfree{iTrl}(jCh,:)]))];
        end
    end
end

switch cfg.criterion.center
    case 'median'
        centerfun = @(x) nanmedian(x,2);
    case 'mean'
        centerfun = @(x) nanmean(x,2);
end
cfg.thres.center = cell2mat(cellfun(@(x) centerfun(x), vardata,'UniformOutput',0));

if strcmp(cfg.criterion.var,'centerMAD')
    varfun = @(x) mad(x,1,2,'omitnan');
elseif any(contains(cfg.criterion.var,{'centerSD','scaleEnvSD','scaleFiltSD'}))
    varfun = @(x) nanstd(x,1,2);
else
    varfun = @(x) [];
end
cfg.thres.variance = cell2mat(cellfun(@(x) varfun(x),vardata,'UniformOutput',0));

if any(contains(cfg.criterion.var,{'centerMAD','centerSD'}))
    thresfun = @(x,y,z) x + (y .* z);
elseif any(contains(cfg.criterion.var,{'scaleEnvSD','scaleFiltSD'}))
    thresfun = @(x,y,z) y .* z;
elseif strcmp(cfg.criterion.var,'scaleCenter')
    thresfun = @(x,y,z) x .* z;
end

if strcmp(cfg.criterion.var,'percentile')
    cfg.thres.main          = cellfun(@(x) prctile(x, cfg.criterion.scale,2), vardata);
    cfg.thres.upperCutoff   = cellfun(@(x) prctile(x, cfg.paramOpt.upperCutoff,2), vardata);

    if cfg.paramOpt.scndThres > 0
        cfg.thres.second = cellfun(@(x) prctile(x, cfg.paramOpt.scndThres,2), vardata);
    end
else
    cfg.thres.main           = thresfun(cfg.thres.center, cfg.thres.variance, cfg.criterion.scale);
    cfg.thres.upperCutoff    = thresfun(cfg.thres.center, cfg.thres.variance, cfg.paramOpt.upperCutoff);

    if cfg.paramOpt.scndThres > 0
        cfg.thres.second = thresfun(cfg.thres.center, cfg.thres.variance, cfg.paramOpt.scndThres);
    end
end
    
clear tmpThresEnv


%% If upperCutoff < Inf set samples upperCutoff to Nan and re-calculate the threshold
if cfg.paramOpt.upperCutoff < Inf
    for iCh = 1 : numThresCh
        vardata{iCh}(vardata{iCh} > cfg.thres.upperCutoff(iCh)) = nan;
    end
    
    cfg.thres.center      = cell2mat(cellfun(@(x) centerfun(x), vardata,'UniformOutput',0));
    cfg.thres.variance    = cell2mat(cellfun(@(x) varfun(x),vardata,'UniformOutput',0));
    
    if strcmp(cfg.criterion.var,'percentile')
        cfg.thres.main = cellfun(@(x) prctile(x, cfg.criterion.scale,2), vardata);
        
        if cfg.paramOpt.scndThres > 0
            cfg.thres.second = cellfun(@(x) prctile(x, cfg.paramOpt.scndThres,2), vardata);
        end
    else
        cfg.thres.main = thresfun(cfg.thres.center, cfg.thres.variance, cfg.criterion.scale);
        
        if cfg.paramOpt.scndThres > 0
            cfg.thres.second = thresfun(cfg.thres.center, cfg.thres.variance, cfg.paramOpt.scndThres);
        end
    end
end    


%% If desired calculate average threshold
if strcmp(cfg.param.thresType, 'average')
    meanThres       = mean(cfg.thres.main,1);
    cfg.thres.main  = repmat(meanThres,numThresCh,1);
    
    if cfg.paramOpt.scndThres > 0
        meanThres           = mean(cfg.thres.second,1);
        cfg.thres.second    = repmat(meanThres,numThresCh,1);
    end
    
    clear meanThres
end


%% Prepare new artifact filter for channels specified for the detection
bnryArtfree = cellfun(@(x) ones(numDetectCh,size(x,2),'logical'),inData.trial,'UniformOutput',0);

if isfield(inData, 'artifacts')
    for iTrl = 1 : numTrl
        for jCh = 1 : numDetectCh
            currCh = ismember(inData.label,cfg.detectCh{jCh});
            
            % @ LDK
            if isempty(inData.artifacts{currCh,iTrl})
                continue
            end
            
            tmpArt = inData.artifacts{currCh,iTrl} + artfctPad;
            tmpArt(tmpArt(:,1) < 1,1)                      = 1;                         %% Ensure padding does not...
            tmpArt(tmpArt(:,2) > inData.sampleinfo(1,2),2) = inData.sampleinfo(1,2);    %% exceed data range
            
            for iArt = 1 : size(tmpArt,1)
                bnryArtfree{iTrl}(jCh,tmpArt(iArt,1):tmpArt(iArt,2)) = 0;
            end
            
            clear tmpArt
        end
    end
end


%% Prepare data to detect events
tmpDetect = [];

for iTrl = 1 : numTrl                                       %% loop over trials
    for jCh = 1 : numDetectCh                               %% loop over channels
        if any(ismember(cfg.thresCh, cfg.detectCh{jCh}))    %% re-use previous filtered data if possible
            tmpDetect.trial{iTrl}(jCh,:) = tmpThres{iTrl}(ismember(cfg.thresCh, cfg.detectCh{jCh}),:);
        else
            fprintf('Channel %d/%d (%s)\n', jCh,numDetectCh, cfg.detectCh{jCh});
            currCh = ismember(inData.label, cfg.detectCh{jCh});
            switch cfg.filtType
                case 'fir'
                    tmpDetect.trial{iTrl}(jCh,:) = ft_preproc_bandpassfilter(inData.trial{iTrl}(currCh,:), fsample, cfg.param.bpfreq, 3*fix(fsample/cfg.param.bpfreq(1))+1, 'fir', 'twopass');
                case 'but'
                    tmpHP = ft_preproc_highpassfilter(inData.trial{iTrl}(currCh,:), fsample, cfg.param.bpfreq(1), 5, 'but', 'twopass','reduce');
                    tmpDetect.trial{iTrl}(jCh,:) = ft_preproc_lowpassfilter(tmpHP, fsample, cfg.param.bpfreq(2), 5, 'but', 'twopass','reduce');
                    
                    clear tmpHP
            end
        end
    end
end

clear tmpThres

%% calculate envelope of bandpass filtered signal
tmpDetectEnv = [];
switch cfg.param.envType
    case 'rms'
        tmpDetectEnv.trial = cellfun(@(x) envelope(x',round(cfg.param.envWin * fsample),'rms')',tmpDetect.trial,'UniformOutput',0);
    case 'env'
        tmpDetectEnv.trial = cellfun(@(x) envelope(x',round(cfg.param.envWin * fsample),'analytic')',tmpDetect.trial,'UniformOutput',0);
end

if cfg.paramOpt.smoothEnv > 0       %% Smooth envelope if requested
    if ~verLessThan('matlab','9.2')
        tmpDetectEnv.trial = cellfun(@(x) smoothdata(x,2,'movmean',round(cfg.paramOpt.smoothEnv * fsample)),tmpDetectEnv.trial,'UniformOutput',0);
    else
        for iTrl = 1 : numTrl
            for jCh = 1 : numDetectCh
                tmpDetectEnv.trial{iTrl}(jCh,:) = smooth(tmpDetectEnv.trial{iTrl}(jCh,:),round(cfg.paramOpt.smoothEnv * fsample));
            end
        end
    end
end


%% detect crossings
supThres = cell(1,numTrl);
for iTrl = 1 : numTrl % loop over trials
    for jCh = 1 : numDetectCh
        currCh = ismember(cfg.thresCh, cfg.detectCh{jCh});  % match current detection channel to thresCh vector
        
        supThres{iTrl}(jCh,:) = all([tmpDetectEnv.trial{iTrl}(jCh,:) >= cfg.thres.main(currCh);...
                                     tmpDetectEnv.trial{iTrl}(jCh,:) <= cfg.thres.upperCutoff(currCh);...
                                     bnryStage{iTrl};...
                                     bnryArtfree{iTrl}(jCh,:)]);
    end
end


%% check all event requirements and calculate metrics
for iTrl = 1 : numTrl               % loop over trials
    for jCh = 1 : numDetectCh       % loop over detect channels
        
        cfg.evtIndiv(jCh,iTrl).label   = cfg.detectCh{jCh};                                      % Channel name
        cfg.evtIndiv(jCh,iTrl).tss     = sum(all([bnryStage{iTrl};...                            % time spend asleep (in min), based on artifact-free sleep
                                                  bnryArtfree{iTrl}(jCh,:)]),2) / (fsample * 60);
                
        %% Optional: Second thresholding
        if cfg.paramOpt.scndThres > 0
            dsig    = diff([0 supThres{iTrl}(jCh,:) 0]);
            staIdx  = find(dsig > 0);
            endIdx  = find(dsig < 0) - 1;
            
            currCh = ismember(cfg.thresCh, cfg.detectCh{jCh});
            for kEvt = 1 : size(staIdx,2)
                if max(tmpDetectEnv.trial{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt))) < cfg.thres.second(currCh)
                    supThres{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt)) = 0;
                end
            end
        end
        
        %% Discard intervals not fulfilling minimal length
        dsig    = diff([0 supThres{iTrl}(jCh,:) 0]);
        staIdx  = find(dsig > 0);
        endIdx  = find(dsig < 0)-1;
        tmpLen  = (endIdx-staIdx+1) < round(cfg.criterion.len(1) * fsample);
        
        rmvIdx                      = cell2mat(reshape(arrayfun(@(x,y) x:y, staIdx(tmpLen), endIdx(tmpLen),'UniformOutput',0),1,sum(tmpLen)));
        supThres{iTrl}(jCh,rmvIdx)  = 0;
        
%         for kEvt = 1 : size(staIdx,2)
%             if duration(kEvt) < round(cfg.criterionLen(1,1) * fsample)
%                 supThres{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt)) = 0;
%             end
%         end
        
        %% Optional: Merge intervals closer than specified margin
        if cfg.paramOpt.mergeEvts > 0
            dsig    = diff([0 supThres{iTrl}(jCh,:) 0]);
            staIdx  = find(dsig > 0);
            endIdx  = find(dsig < 0)-1;
            
            for kEvt = 1 : size(endIdx,2)-1 %% Check additionally if gap is without artifacts
                if staIdx(kEvt+1) - endIdx(kEvt) <= round(cfg.paramOpt.mergeEvts * fsample) &&...
                   all(bnryArtfree{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt+1)))
                    supThres{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt+1)) = 1;
                end
            end
        end
        
        %% Optional: Discard events not fulfilling minimal number of cycles
        if cfg.paramOpt.minCycles > 0
            dsig    = diff([0 supThres{iTrl}(jCh,:) 0]);
            staIdx  = find(dsig > 0);
            endIdx  = find(dsig < 0)-1;
            
            switch cfg.paramOpt.minCycles
                case 1
                    currCh  = ismember(inData.label, cfg.detectCh{jCh});
                    currSig = smoothdata(inData.trial{iTrl}(currCh,:),'movmedian',3);
                case 2
                    currSig = smoothdata(tmpDetect.trial{iTrl}(jCh,:),'movmedian',3);
            end
            
            for kEvt = 1 : size(staIdx,2)
                [~,maxidx] = findpeaks(currSig(1,staIdx(kEvt):endIdx(kEvt)));
                [~,minidx] = findpeaks((-1) * currSig(1,staIdx(kEvt):endIdx(kEvt)));
                
                if (numel(maxidx) < cfg.paramOpt.minCyclesNum) || (numel(minidx) < cfg.paramOpt.minCyclesNum)
                    supThres{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt)) = 0;
                end
            end
        end
        
        %% Last iteration to discard spindles longer than specified duration
        numEvt = 0; % initialize event counter
        
        dsig    = diff([0 supThres{iTrl}(jCh,:) 0]);
        staIdx  = find(dsig > 0);
        endIdx  = find(dsig < 0)-1;
        tmpLen  = endIdx-staIdx+1;
        
        for kEvt = 1: size(staIdx,2)    % loop over potential spindles
            if tmpLen(kEvt) <= round(cfg.criterion.len(2)*fsample) && ... % inclusion criterion fullfilled w/o artifacts
               all(all([bnryStage{iTrl}(staIdx(kEvt) - criterionPad(1) : endIdx(kEvt) + criterionPad(2)); ...
                        bnryArtfree{iTrl}(jCh,staIdx(kEvt) - criterionPad(1) : endIdx(kEvt) + criterionPad(2))]))
                
                numEvt = numEvt + 1;
                
                cfg.evtIndiv(jCh,iTrl).staTime(numEvt) = staIdx(kEvt);                                 % event start (in datapoints)
                cfg.evtIndiv(jCh,iTrl).midTime(numEvt) = round(mean([staIdx(kEvt) endIdx(kEvt)]));   % event mid (in datapoints)
                cfg.evtIndiv(jCh,iTrl).endTime(numEvt) = endIdx(kEvt);                                   % event end (in datapoints)
                
                cfg.evtIndiv(jCh,iTrl).stage(numEvt)   = inData.staging{iTrl}(1,staIdx(kEvt));
                
                cfg.evtIndiv(jCh,iTrl).duration(numEvt) = (endIdx(kEvt)-staIdx(kEvt)) / fsample;  % event duration (in seconds)
                
                tmpWin = tmpDetect.trial{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt));
                [minAmp,minIdx] = min(tmpWin);
                [maxAmp,maxIdx] = max(tmpWin);
                cfg.evtIndiv(jCh,iTrl).maxTime(numEvt) = staIdx(kEvt) + maxIdx - 1; % time of event peak (in datapoints)
                cfg.evtIndiv(jCh,iTrl).minTime(numEvt) = staIdx(kEvt) + minIdx - 1; % time of event trough (in datapoints)
                cfg.evtIndiv(jCh,iTrl).minAmp(numEvt) = minAmp;
                cfg.evtIndiv(jCh,iTrl).maxAmp(numEvt) = maxAmp;
                
                tmpWin                                      = tmpDetectEnv.trial{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt));
                [envMaxAmp,envMaxIdx]                       = max(tmpWin);
                cfg.evtIndiv(jCh,iTrl).envMaxAmp(numEvt)    = envMaxAmp;                     % RMS max
                cfg.evtIndiv(jCh,iTrl).envMaxTime(numEvt)   = staIdx(kEvt) + envMaxIdx - 1;  % time of RMS max (in datapoints)
                cfg.evtIndiv(jCh,iTrl).envMean(numEvt)      = mean(tmpWin,2);
                cfg.evtIndiv(jCh,iTrl).envSum(numEvt)       = sum(tmpWin,2);
                
                % Add peaks and troughs
                tmpWin = tmpDetect.trial{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt));
                [~, EvtPeaks]   = findpeaks(tmpWin,staIdx(kEvt):endIdx(kEvt));
                [~, EvtTroughs] = findpeaks((-1) * tmpWin,staIdx(kEvt):endIdx(kEvt));
                                                
                cfg.evtIndiv(jCh,iTrl).peaks{numEvt}    = EvtPeaks;
                cfg.evtIndiv(jCh,iTrl).troughs{numEvt}  = EvtTroughs;
                cfg.evtIndiv(jCh,iTrl).freq(numEvt)     = mean([diff(EvtPeaks),diff(EvtTroughs)]);
                                
            else % duration criterion NOT fullfilled
                supThres{iTrl}(jCh,staIdx(kEvt):endIdx(kEvt)) = 0; % remove unfit events from suprathresh-data
            end
        end
        
        cfg.evtIndiv(jCh,iTrl).numEvt = numEvt;   % Save number of detected events
    end
end

clear tmpDetectEnv tmpDetect supThres


%% Optional: False-positive rejection based on frequency profile
if cfg.doFalseposRjct
    fprintf('----- Event rejection by frequency profile\n');
    cfg.falseposRjct.rejects = cell(numDetectCh,numTrl);
    
    for iTrl = 1 : numTrl
        for jCh = 1 : numDetectCh
            fprintf('Channel %d/%d (%s): ', jCh, numDetectCh,cfg.detectCh{jCh});
            
            tmpTic = tic;
            currCh = ismember(inData.label,cfg.detectCh{jCh});
            
            %--- Segment data
            tfg         = [];
            tfg.trl     = [cfg.evtIndiv(jCh,iTrl).envMaxTime' + round(cfg.falseposRjct.timePad(1) * fsample),...
                           cfg.evtIndiv(jCh,iTrl).envMaxTime' + round(cfg.falseposRjct.timePad(2) * fsample),...
                           ones(cfg.evtIndiv(jCh,iTrl).numEvt,1) * round(cfg.falseposRjct.timePad(1) * fsample)];
            tmpTrls   = ft_redefinetrial(tfg,inData);
            
            %--- Calculate time frequency representation
            tfg             = [];
            tfg.channel     = tmpTrls.label(currCh);
            tfg.taper       = 'hanning';
            tfg.method      = 'mtmconvol';
            tfg.pad         = 'nextpow2';
            tfg.output      = 'pow';
            tfg.keeptrials  = 'yes';
            tfg.foi         = cfg.falseposRjct.tfrFreq;
            tfg.toi         = cfg.falseposRjct.tfrTime;
            tfg.t_ftimwin   = cfg.falseposRjct.tfrWin;
            
            tmpTFR  = ft_freqanalysis(tfg,tmpTrls);
            tmpPow  = squeeze(tmpTFR.powspctrm);                                        %% Note: rpt x freq x time
            tmpTime = arrayfun(@(x) nearest(tmpTFR.time,x),cfg.falseposRjct.avgWin);
            tmpFreq = arrayfun(@(x) nearest(tmpTFR.freq,x),cfg.falseposRjct.freqlim);
            
            %--- Perform event rejection
            cfg.falseposRjct.rejects{jCh,iTrl} = ones(size(tmpPow,1),1,'logical');
            
            tmpPow = squeeze(sum(tmpPow(:,:,tmpTime(1):tmpTime(2)),3));
            tmpMax = max(tmpPow,[],2);                                      %% Determine maximum per trial
            tmpPow = tmpPow ./ repmat(tmpMax,1,size(tmpPow,2));             %% Normalise by maximum value
  
            for iEvt = 1 : size(tmpPow,1)
                [~, tmpPks,~,tmpProm] = findpeaks(tmpPow(iEvt,:));
                
                hazMax = find(tmpPks >= tmpFreq(1) & tmpPks <= tmpFreq(2));
                if numel(hazMax) > 0 && any(tmpProm(hazMax) > cfg.falseposRjct.prominence)
                    cfg.falseposRjct.rejects{jCh,iTrl}(iEvt) = 0;
                end
            end
            
            fprintf(' reject %d of %d (%.2f) - took %.2f s\n',...
                     sum(cfg.falseposRjct.rejects{jCh,iTrl}),...
                     size(tmpPow,1),...
                     100 * sum(cfg.falseposRjct.rejects{jCh,iTrl}) / size(tmpPow,1),...
                     toc(tmpTic));
            
            clear tmpTrls tmpTFR tmpPow

        end
    end
end


%% add summary statistics to cfg.evtSummary
for iTrl = 1 : numTrl % loop over trials
    for jCh = 1 : numDetectCh
        cfg.evtSummary(jCh,iTrl).label    = cfg.detectCh{jCh};
        cfg.evtSummary(jCh,iTrl).numEvt   = cfg.evtIndiv(jCh,iTrl).numEvt;
        cfg.evtSummary(jCh,iTrl).tss      = cfg.evtIndiv(jCh,iTrl).tss;
        
        if cfg.evtSummary(jCh,iTrl).numEvt > 0
            if cfg.doFalseposRjct
                tmpIdx = ~cfg.falseposRjct.rejects{jCh,iTrl};
            else
                tmpIdx = ones(cfg.evtSummary(jCh,iTrl).numEvt,1,'logical');
            end
            
            cfg.evtSummary(jCh,iTrl).meanDuration   = mean(cfg.evtIndiv(jCh,iTrl).duration(tmpIdx),2);
            cfg.evtSummary(jCh,iTrl).meanFreq       = mean(cfg.evtIndiv(jCh,iTrl).freq(tmpIdx),2);
            cfg.evtSummary(jCh,iTrl).meanMin        = mean(cfg.evtIndiv(jCh,iTrl).minAmp(tmpIdx),2);
            cfg.evtSummary(jCh,iTrl).meanMax        = mean(cfg.evtIndiv(jCh,iTrl).maxAmp(tmpIdx),2);
            cfg.evtSummary(jCh,iTrl).meanEnvMax     = mean(cfg.evtIndiv(jCh,iTrl).envMaxAmp(tmpIdx),2);
            cfg.evtSummary(jCh,iTrl).meanEnvMean    = mean(cfg.evtIndiv(jCh,iTrl).envMean(tmpIdx),2);
            cfg.evtSummary(jCh,iTrl).meanEnvSum     = mean(cfg.evtIndiv(jCh,iTrl).envSum(tmpIdx),2);
        end
    end
end
end