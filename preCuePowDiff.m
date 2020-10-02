%% QUESTIONS
%  Should I orthogonalize?
%  If so, should I normalize before or after orthogonalization?
%  I am now comparing only idx and ndx power for bundles with indexing activity
%  Should I mean over all channel? I think meaning should have an effect
%  although I orthogonalize. However, there is an argument that I should no
%  longer need to mean over all channels since they do not share variance
%  anymore

function preCuePowDiff % does not consider artefacts yet
lfpDir = dir('X:\Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_noSPKINT.mat');  % CHANGED TO NO SPK INT
% artFold = 'Z:\hanslmas-ieeg-compute\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS
artFold = 'X:\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS
load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')
counter = 1;
skippedDat = [];
hz = linspace(0, 1000, 1001);
ortho = 1;

ndxPow    = [];
idxPow    = [];
% diffPow   = [];
% bundleVar = [];
for spk = 1 : length(allSpks)
    
    %% SKIP REPEATING MICROWIRE BUNDLES
    if spk > 1
        if and(and(contains(allSpks(spk-1).bidsID, allSpks(spk).bidsID), strcmp(allSpks(spk-1).sesh, allSpks(spk).sesh) ), strcmp(allSpks(spk).bundlename, allSpks(spk-1).bundlename)); % same subject + session
            continue
        end
    end
    
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
    curBund = allSpks(spk).bundlename;
    allBund = {allSpks.bundlename};
    
    %% LOAD IN THE LFP-DATA
    try
        load([lfpDir(1).folder, filesep, bidsID, '_', sesh, '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    catch
        skippedDat = [skippedDat; {bidsID} {sesh}];
        disp('skipDat');
        error('No files should be skipped anymore.')
        continue
    end
    
    %% normalize LFP variance
    microLFP.trial = (microLFP.trial - mean(microLFP.trial,2)) / std(microLFP.trial,0,2);
    
    %% ORTHOGONALIZE LFP PER BUNDLE
    switch ortho
        case 0
            dataBundLab = cellfun(@(x) x(1:end-1) , data.label, 'un', 0); % extract the bundle of interest
        case 1
            data = orthogonVec(data); % with orthogonalization
            dataBundLab = cellfun(@(x) x(1:end-8) , data.label, 'un', 0); % during orthogonVec the labels are extended so it takes more to get the bundle definition
    end
    
    loadBunds    = strcmp(curBund, dataBundLab);
    
    cfg          = [];
    cfg.channel  = data.label(loadBunds);     % all 8 channels
    microLFP     = ft_selectdata(cfg, data);  % select
    
    
    
    %% WHICH TRIALS IN THAT BUNDLE ARE INDEXED?
    sameBund = and(and(contains({allSpks.bidsID}, bidsID), strcmp({allSpks.sesh}, sesh) ), contains(allBund, curBund)); % same subject + session + bundle
    idxTrl   = any(vertcat(allSpks(sameBund).idxTrl),1); % these are the trials that are indexed in that wire
    
    encTrig  = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1])*1000);
    
    %% REDEFINE TRIALS ACCORDING TO encTrig
    cfg          = [ ];
    cfg.trl      = [encTrig-1000 encTrig zeros(size(encTrig))];
    microLFPtrl  = ft_redefinetrial(cfg, microLFP);
    
    curNdx_allChan  = [];
    curIdx_allChan  = [];
    curIdx_chanMean = [];
    curNdx_chanMean = [];
    for trl = 1 : length(encTrig)
        %         trlLFP = microLFP(:,encTrig(trl)-1000:encTrig(trl)); % before FT
        %         trlPow = abs(fft(trlLFP)); % before FT
        %         trlPow = mean(trlPow, 1);
        
        % FIELDTRIP POWER
        cfg = [];
        cfg.output    = 'pow';
        cfg.channel   = 'all';
        cfg.method    = 'mtmfft';
        cfg.pad       = 4; % padding to increase frequency resolution
        cfg.foi       = [1:0.25:200];
        cfg.taper     = 'hanning';
        cfg.trials    = trl;
        trlPow        = ft_freqanalysis(cfg, microLFPtrl);
        
        switch idxTrl(trl)
            case 0 % trial IS NOT being indexed
                ndxPow          = [ndxPow;                  trlPow.powspctrm];
                curNdx_chanMean = [curNdx_chanMean; mean(trlPow.powspctrm,1)]; % SHOULD I MEAN OVER CHANNELS HERE?
                curNdx_allChan  = [curNdx_allChan,        trlPow.powspctrm,1];
            
            case 1 % trial IS being indexed
                idxPow          = [idxPow;                  trlPow.powspctrm];
                curIdx_chanMean = [curIdx_chanMean; mean(trlPow.powspctrm,1)]; % SHOULD I MEAN OVER CHANNELS HERE?
                curIdx_allChan  = [curIdx_allChan,        trlPow.powspctrm,1];

%                 bundleVar = [bundleVar {trlPow.powspctrm}]; % I want to see if there is a difference within a bundle
        end
        
    end
    
    if ~isempty(curIdx)
        %         diffPow   = [diffPow; mean(curIdx,1)-mean(curNdx,1)];
        
        %% THE MEAN OVER ALL 8 MICROWIRES
        bundlePow_chanMean.idx = [ bundlePow_chanMean.idx {curIdx_chanMean} ];
        bundlePow_chanMean.ndx = [ bundlePow_chanMean.ndx {curNdx_chanMean} ];
        
        %% ALL 8 MICROWIRES AS RANDOM EFFECT
        bundlePow_allChan.idx  = [ bundlePow_allChan.idx  {curIdx_allChan}  ];
        bundlePow_allChan.ndx  = [ bundlePow_allChan.ndx  {curNdx_allChan}  ];
        
    end
    
end

hz = trlPow.freq;
save('X:\Luca\data\allSbj\preCuePowDiff_orthNorm.mat', 'idxPow', 'ndxPow', 'bundlePow', 'hz')

%% LOOKING AT POWER AS A RANDOM EFFECT NOW (CAN ALSO USE diffPow AS A WITHIN BUNDLE POWER DIFFERENCE!

%% FREQUENCY BANDS
freq.delta  = hz<4;
freq.theta  = hz>= 4 & hz< 8;
freq.alphaL = hz>= 8 & hz<12;
freq.alphaH = hz>=12 & hz<16; 
freq.beta   = hz>=16 & hz<30;

%% EMPIRICAL POWER DIFFERENCE
powDiff = median(idxPow,1)- median(ndxPow,1);
powDiff = [sum(powDiff(:,delta),2) sum(powDiff(:,theta),2) sum(powDiff(:,alphaL),2) sum(powDiff(:,alphaH),2) sum(powDiff(:,beta),2)];   

%% PERMUTATIOM
allPow = [idxPow; ndxPow];
numIdx = size(idxPow, 1); % number of indexed trials
nperm  = 10000;

powDiffPerm = zeros(nperm, length(hz));
for perm = 1:nperm
    disp(perm)
    randDx     = randperm(size(allPow, 1)); % shuffle power spectra
    randIdx    = randDx(1:numIdx); % permuted idx power spectra
    randNdx    = randDx(numIdx+1:end); % permuted ndx power spectra
    
    idxPowPerm = allPow( randIdx, :);
    ndxPowPerm = allPow( randNdx, :);
    
    powDiffPerm(perm,:) = median(idxPowPerm,1) - median(ndxPowPerm,1); 
end

% BINNING
save('Z:\hanslmas-ieeg-compute\Luca\data\allSbj\powDiffPerm.mat', 'powDiffPerm');
powDiffPermBin = [sum(powDiffPerm(:,delta),2) sum(powDiffPerm(:,theta),2) sum(powDiffPerm(:,alphaL),2) sum(powDiffPerm(:,alphaH),2) sum(powDiffPerm(:,beta),2)];   

%% PLOTTING
figure(1); clf; hold on;
% plot([0 5], [0 0], 'linew', 2, 'color', 'k')

plot(powDiff, 'linew', 2, 'color', 'r');
plot(prctile(powDiffPermBin, 100-(2.5/5)), 'linew', 2, 'color', 'k');
plot(prctile(powDiffPermBin, (2.5/5)), 'linew', 2, 'color', 'k');

xticks([1:5])
xticklabels({'Delta', 'Theta', 'Low Alpha', 'High Alpha', 'Beta'})
ylabel('Power Indexed Trials Minus Non-Indexed Trials')

scatter(1, -20, 100, 'k', 'filled', 'p')
% scatter(2, -20, 100, 'k', 'filled', 'p')
% scatter(3, -20, 100, 'k', 'filled', 'p')
% scatter(4, -20, 100, 'k', 'filled', 'p')
% scatter(5, -20, 100, 'k', 'filled', 'p')

legend({'Empirical Power Difference', '97.5- and 2.5-Percentile'})
title({'Power Difference Between Indexed And Non-Indexed Trials ','(Pre-Cue & Corrected)'})


% %% LOOK AT WITHIN BUNDLE VARIANCE?
% figure(1)
% for idx = 1:length(bundleVar)
%     clf
%     for mw = 1:8
%         subplot(2,4,mw)
%         plot(hz, bundleVar{idx}(mw,:), 'b', 'linew', 1)
%         % xlim([0 100])
%         drawnow
%         scale(mw,:) = get(gca, 'YLim');
%     end
%     
%     scale = max(scale, [], 1);
%     for mw = 1:8
%         subplot(2,4,mw)
%         set(gca, 'YLim', [0 scale(2)]);
%     end
%     
%     pause(5)
% end

end % END OF FUNCTION
