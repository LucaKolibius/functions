function TFtoi3(trigCode,TW)
whereAmI(0);
global prePath;

savenam = ['TFtoi3_trigCode_', num2str(trigCode)];

lfpDir = dir([prePath, 'Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_*.mat']);  % spkInt in laptop and noSpkInt on local
load([prePath, 'Luca\data\allSbj\allSpksHZ.mat'], 'allSpks')

allSUPowLow.idx   = [];
allSUPowLow.ndx   = [];
allSUPowHigh.idx  = [];
allSUPowHigh.ndx  = [];


prevLFP = [];
for spk = 1 : length(allSpks)
%     try
        if any(isnan(allSpks(spk).idxTrlSingLw)) & any(isnan(allSpks(spk).idxTrlSingHi))
            continue
        end
        
        %% GET: bidsID + sesh
        bidsID    = allSpks(spk).bidsID;
        sesh      = allSpks(spk).sesh;
        favChanHi = allSpks(spk).favChanHigh;
        favChanLo = allSpks(spk).favChanLow;
        idxTrlLw  = allSpks(spk).idxTrlSingLw;
        idxTrlHi  = allSpks(spk).idxTrlSingHi;
        encTrig   = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx, trigCode)*1000);
        bund      = allSpks(spk).bundlename;
        
        %% LOAD IN THE LFP-DATA
        newLFP = [lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'];
        if ~strcmp(newLFP, prevLFP)
            clear data
            load(newLFP, 'data');
            prevLFP = newLFP;
        end
        
        %% SELECT BUNDLE THAT REFLECTS SU INPUT
        allBund = cellfun(@(x) x(1:end-1), data.label, 'un', 0);
        bundIdx = strcmp(allBund, bund);
        
        cfg         = [];
        cfg.channel = data.label(bundIdx);
        selData  = ft_selectdata(cfg, data);
        
        %% REDEFINE TRIALS ACCORDING TO encTrig (WITH WINGS FOR TF)
        switch size(encTrig,2)
            case 1 % snippet in relation to one timepoint (cuelocked or resplocked)
                cfg          = [ ];
                cfg.trl      = [encTrig-TW(1)-3000 encTrig-1+TW(2)+3000 zeros(size(encTrig))];
                microLFPex   = ft_redefinetrial(cfg, selData);
                
                % WITHOUT WINGS (FOR AR)
                cfg          = [ ];
                cfg.trl      = [encTrig-TW(1) encTrig-1+TW(2) zeros(size(encTrig))];
                microLFP     = ft_redefinetrial(cfg, selData);
                
            case 2 % snippet in relation to two timepoints (cuelocked to resplocked / whole trial)
                cfg          = [ ];
                cfg.trl      = [encTrig(:,1)-TW(1)-3000 encTrig(:,2)-1+TW(2)+3000 zeros(size(encTrig,1))];
                microLFPex   = ft_redefinetrial(cfg, selData);
                
                % WITHOUT WINGS (FOR AR)
                cfg          = [ ];
                cfg.trl      = [encTrig(:,1)-TW(1) encTrig(:,2)-1+TW(2) zeros(size(encTrig,1))];
                microLFP     = ft_redefinetrial(cfg, selData);
        end
        
        
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
        allPowLow   = [];
        allPowHigh  = [];
        
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
        tempVar(spk).allPowLow = allPowLow;
        %     allSpks(spk).allPowLow = allPowLow;
        
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
        %     allSpks(spk).allPowHigh = allPowHigh;
        tempVar(spk).allPowHigh = allPowHigh;
        
        
        %% UPDATE PROGRESS
        lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
        lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
        lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
        lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
        lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
        
%     catch
%         %         save(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\', savenam, '_partial.mat'], 'allSpks', '-v7.3')
%         save(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\TFtoi\', savenam,  'spk_' num2str(spk), '_partial.mat'], 'allSpks', '-v7.3')
%         
%     end
end % END OF SU LOOP


% freqRes = trlPow.freq;
% save(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\', savenam, '.mat'], 'allSpks', '-v7.3')
save(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\TFtoi\', savenam,  'spk_' num2str(spk), '.mat'], 'allSpks', '-v7.3')

end % END OF FUNCTION