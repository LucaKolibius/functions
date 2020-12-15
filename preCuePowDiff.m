function preCuePowDiff
global prePath;

lfpDir = dir([prePath, 'Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_*.mat']);  % spkInt in laptop and noSpkInt on local
load([prePath, 'Luca\data\allSbj\allSpksHZ.mat'], 'allSpks')

skippedDat = [];
mtrack = [];

allSUPow.idx  = [];
allSUPow.ndx  = [];

for spk = 1 : length(allSpks)
    
    if any(isnan(allSpks(spk).idxTrlSing))
        continue
    end
    
    if isempty(allSpks(spk).favChan)
        error('spPhase_encRet fully completed!')
    end
    
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
    favChan = allSpks(spk).favChan;
    idxTrl  = allSpks(spk).idxTrlSing;     % WHICH TRIALS DOES THAT SU INDEX?
    encTrig = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1])*1000);
    bund = allSpks(spk).bundlename;
    
    %% LOAD IN THE LFP-DATA
    load([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    
    %% SELECT BUNDLE THAT REFLECTS SU INPUT (see spkPhase_encRet)
    allBund = cellfun(@(x) x(1:end-1), data.label, 'un', 0);
    bundIdx = strcmp(allBund, bund);
    
    cfg         = [];
    cfg.channel = data.label(bundIdx);
    microLFP    = ft_selectdata(cfg, data);
    
    
    %% REDEFINE TRIALS ACCORDING TO encTrig
    cfg          = [ ];
    cfg.trl      = [encTrig-1000-10000-500 encTrig-1+10000 zeros(size(encTrig))]; % 500 is the BL, rest are wings
    microLFP     = ft_redefinetrial(cfg, microLFP);
    
    
    %% CONSIDER BANDPASS-FILTERING 1-300HZ!!
    
    %% ONLY DEMEAN & ORTHOGONALIZE
    % normalize LFP variance
    %     data.trial = {(data.trial{1} - mean(data.trial{1},2)) ./ std(data.trial{1},0,2)};
    
    cfg        = [];
    cfg.demean = 'yes';
    microLFP   = ft_preprocessing(cfg, microLFP);
    microLFP   = orthogonVec(microLFP);           % orthogonalize
    
    
    
    
    idxTrlPow = [];
    ndxTrlPow = [];
    for trl = 1 : length(encTrig)
        
        cfg = [];
        cfg.output    = 'pow';
        cfg.channel   = 'all';
        cfg.foi       = logspace(log10(1),log10(140));
        cfg.trials    = trl;
        cfg.toi       = 'all';
        cfg.width     = 6; %12; % 12 for HF (>40hz)
        cfg.method    = 'wavelet';
        
        trlPow        = ft_freqanalysis(cfg, microLFP);
        freqRes       = trlPow.freq;
        trlPow        = trlPow.powspctrm;
        
        %% AR
        [isAR,~]      = iqrAR(microLFP.trial{trl},0);
        isAR          = squeeze(isAR(:,1:length(freqRes),:));
        trlPow(isAR)  = nan;
        
        %% Chose Input Channel per Frequency
        temp = [];
        for freq = 1:length(freqRes)
            realFreq       = freqRes(freq);  % frequency of interest
            freqDiff       = abs(realFreq-cfg.foi); % find difference from frequency of interest and frequencies for which I have favChan
            [~, bestMatch] = min(freqDiff); % find index of frequency with minimum difference
            temp(freq,:)   = trlPow(favChan(bestMatch),freq,:);
        end
        trlPow             = temp;
        
        %% Cut Wings
        trlPow(:,1:10000)      = [];
        trlPow(:,end-9999:end) = [];
%         imagesc(trlPow)

           %% NORMALIZE TF PLOT
    baseline = trlPow(:,1:500);
    trlPow   = trlPow(:,501:end);
    
    meanBL = nanmean(baseline,2);
    stdBL  = nanstd(baseline,0,2);
    trlPow = (trlPow-meanBL)./stdBL;
        
        switch idxTrl(trl)
            case 0
                ndxTrlPow = [ndxTrlPow; {trlPow}];
            case 1
                idxTrlPow = [idxTrlPow; {trlPow}]; 
        end
        
    end
    
    if sum(idxTrl) > 0
        
        allSUPow.idx  = [ allSUPow.idx  {idxTrlPow}  ];
        allSUPow.ndx  = [ allSUPow.ndx  {ndxTrlPow}  ];
        
        addthis = [spk; length(allSUPow.idx)];
        mtrack  = [mtrack, addthis];
    end
    
    %% UPDATE PROGRESS
    lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).', spk, length(allSpks), spk/length(allSpks)*100);
    %     fprintf(repmat('\b',1,lineLength))
    
end
freqRes = trlPow.freq;
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\preCuePowDiff_orthDeMea.mat', 'allSUPow', 'hz')

end % END OF FUNCTION