%% USE TFtoi3 INSTEAD AND ADAPT THE TRIGGER
function TFtoi3_resp
whereAmI(0);
global prePath;

lfpDir = dir([prePath, 'Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_*.mat']);  % spkInt in laptop and noSpkInt on local
load([prePath, 'Luca\data\allSbj\allSpksHZ.mat'], 'allSpks')

allSUPowLow.idx   = [];
allSUPowLow.ndx   = [];
allSUPowHigh.idx  = [];
allSUPowHigh.ndx  = [];

for spk = 1 : length(allSpks)
    try
    %     if any(isnan(allSpks(spk).idxTrlSing))
    %         continue
    %     end
    
    %% GET: bidsID + sesh
    bidsID    = allSpks(spk).bidsID;
    sesh      = allSpks(spk).sesh;
    favChanHi = allSpks(spk).favChanHigh;
    favChanLo = allSpks(spk).favChanLow;
    idxTrlLw  = allSpks(spk).idxTrlSingLw;
    idxTrlHi  = allSpks(spk).idxTrlSingHi;
    encTrig   = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[3])*1000);
    bund      = allSpks(spk).bundlename;
    
    %% LOAD IN THE LFP-DATA
    load([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    
    %% SELECT BUNDLE THAT REFLECTS SU INPUT
    allBund = cellfun(@(x) x(1:end-1), data.label, 'un', 0);
    bundIdx = strcmp(allBund, bund);
    
    cfg         = [];
    cfg.channel = data.label(bundIdx);
    microLFPex  = ft_selectdata(cfg, data);
    
    
    %% REDEFINE TRIALS ACCORDING TO encTrig (WITH WINGS FOR TF)
    cfg          = [ ];
    cfg.trl      = [encTrig-5000-3000 encTrig-1+5000+3000 zeros(size(encTrig))]; % -1s to +5s plus wings
    microLFPex   = ft_redefinetrial(cfg, microLFPex);
    
    % WITHOUT WINGS (FOR AR)
    cfg          = [ ];
    cfg.trl      = [encTrig-5000 encTrig-1+5000 zeros(size(encTrig))];
    microLFP     = ft_redefinetrial(cfg, microLFPex);
    
    
    %% CONSIDER BANDPASS-FILTERING 1-300HZ!!
    
    %% ONLY DEMEAN & ORTHOGONALIZE
    % normalize LFP variance
    %     data.trial = {(data.trial{1} - mean(data.trial{1},2)) ./ std(data.trial{1},0,2)};
    
    cfg         = [];
    cfg.demean  = 'yes';
    microLFPex  = ft_preprocessing(cfg, microLFPex);
    microLFPex  = orthogonVec(microLFPex);           % orthogonalize
    
    cfg         = [];
    cfg.demean  = 'yes';
    microLFP    = ft_preprocessing(cfg, microLFP);
    microLFP    = orthogonVec(microLFP);           % orthogonalize
    
    %% HIGH AND LOW FREQUENCY INPUT - ONE CHANNEL LFP
    allPowLow  = [];
    allPowHigh = [];
    
    %% LOW
    for trl = 1 : length(encTrig)
        
        if any(isnan(idxTrlLw))
            continue
        end
        
        cfg = [];
        cfg.output    = 'pow';
        cfg.channel   = microLFPex.label(favChanLo);
        %         cfg.foi       = logspace(log10(1),log10(140));
        cfg.foi       = 1:1:45;
        cfg.trials    = trl;
        cfg.toi       = 'all';
        cfg.width     = 6; %12; % 12 for HF (>40hz)
        cfg.method    = 'wavelet';
        
        trlPowLow     = ft_freqanalysis(cfg, microLFPex);
        freqResLow    = trlPowLow.freq;
        trlPowLow     = squeeze(trlPowLow.powspctrm);
       
         %% Cut Wings
        trlPowLow(:,1:3000)        = [];
        trlPowLow(:,end-2999:end)  = [];
        
        %% AR
        [isAR,~]             = iqrAR(microLFP.trial{trl},0);
        isARLow              = squeeze(isAR(favChanLo,1:length(freqResLow),:));
        trlPowLow(isARLow)   = NaN;
        
        %         %% Chose Input Channel per Frequency
        %         temp = [];
        %         for freq = 1:length(freqRes)
        %             realFreq       = freqRes(freq);  % frequency of interest
        %             freqDiff       = abs(realFreq-cfg.foi); % find difference from frequency of interest and frequencies for which I have favChan
        %             [~, bestMatch] = min(freqDiff); % find index of frequency with minimum difference
        %             temp(freq,:)   = trlPow(favChan(bestMatch),freq,:);
        %         end
        %         trlPow             = temp;
        
       
        
        %            %% NORMALIZE TF PLOT USING BASELINE (OLD, STANDARDIZE OVER TRIALS NOW)
        %     baseline = trlPow(:,1:500);
        %     trlPow   = trlPow(:,501:end);
        %
        %     meanBL = nanmean(baseline,2);
        %     stdBL  = nanstd(baseline,0,2);
        %     trlPow = (trlPow-meanBL)./stdBL;
        
        %% FROM TF TO POWERSPCTRM
%         trlPowLow  = mean(trlPowLow,2)';
        
        %% SAVE TRL POWER TO ALL POWER
        allPowLow  = cat(3, allPowLow, trlPowLow);
    end
    allSpks(spk).allPowLow = allPowLow;
    
    %% HIGH
    for trl = 1:length(encTrig)
        
        if any(isnan(idxTrlHi))
            continue
        end
        
        cfg = [];
        cfg.output    = 'pow';
        cfg.channel   = microLFPex.label(favChanHi);
        %         cfg.foi       = logspace(log10(1),log10(140));
        cfg.width    = 12;
        cfg.foi       = 70:1:150;
        cfg.trials    = trl;
        cfg.toi       = 'all';
        cfg.method    = 'wavelet';
        
        trlPowHigh    = ft_freqanalysis(cfg, microLFPex);
        freqResHigh   = trlPowHigh.freq;
        trlPowHigh    = squeeze(trlPowHigh.powspctrm);
        
                %% Cut Wings
        trlPowHigh(:,1:3000)       = [];
        trlPowHigh(:,end-2999:end) = [];

        %% AR
        [isAR,~]             = iqrAR(microLFP.trial{trl},0);
        isARHigh             = squeeze(isAR(favChanHi,1:length(freqResHigh),:));
        trlPowHigh(isARHigh) = NaN;
        
        %         %% Chose Input Channel per Frequency
        %         temp = [];
        %         for freq = 1:length(freqRes)
        %             realFreq       = freqRes(freq);  % frequency of interest
        %             freqDiff       = abs(realFreq-cfg.foi); % find difference from frequency of interest and frequencies for which I have favChan
        %             [~, bestMatch] = min(freqDiff); % find index of frequency with minimum difference
        %             temp(freq,:)   = trlPow(favChan(bestMatch),freq,:);
        %         end
        %         trlPow             = temp;
        
        
        %            %% NORMALIZE TF PLOT USING BASELINE (OLD, STANDARDIZE OVER TRIALS NOW)
        %     baseline = trlPow(:,1:500);
        %     trlPow   = trlPow(:,501:end);
        %
        %     meanBL = nanmean(baseline,2);
        %     stdBL  = nanstd(baseline,0,2);
        %     trlPow = (trlPow-meanBL)./stdBL;
        
       
        %% SAVE TRL POWER TO ALL POWER
        allPowHigh  = cat(3, allPowHigh, trlPowHigh);
        
    end
    allSpks(spk).allPowHigh = allPowHigh;
    
    
    %% UPDATE PROGRESS
lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).', spk, length(allSpks), spk/length(allSpks)*100);
lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).', spk, length(allSpks), spk/length(allSpks)*100);
lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).', spk, length(allSpks), spk/length(allSpks)*100);
lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).', spk, length(allSpks), spk/length(allSpks)*100);
lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).', spk, length(allSpks), spk/length(allSpks)*100);

%     fprintf(repmat('\b',1,lineLength))
    
catch
        save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\TFtoi3_respLocked_partial.mat', 'allSpks', '-v7.3')
    end

end % END OF SU LOOP


% freqRes = trlPow.freq;
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\TFtoi3_respLocked.mat', 'allSpks', '-v7.3')

end % END OF FUNCTION