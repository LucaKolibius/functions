%% FavChan Analyse
% is allSpksHZaftTrans (after transfer from bham) the same as the one I
% have on castles?? I did some calculations until after midnight on the 6th
% early 6th of november, but htis variable is from early morning thursday.
% I also should have a allSpksHZ_phs variable here that is missing

clear
server = 0;
switch server
    case 0     % RUN ON MY OWN LAPTOP
        addpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\toolboxes\fieldtrip-20200603');
        addpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\toolboxes\CircStat2012a');
        ft_defaults
        allMicro = '\\analyse4.psy.gla.ac.uk\project0309\Luca\data\microLFP\';
        load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat', 'allSpks');
        
    case 1         % RUN ON THE SERVER
end

prevLFPname = [];
mkPlots = 0;
for su = 252: size(allSpks,2)
    
    disp(su)
    
    spks = allSpks(su).spks;
    spks = round(spks);
    
    %     if size(spks,1) > 1000
    %         idx = randperm(size(spks, 1))
    %         spks = spks(idx);
    %         spks = spks(1:1000);
    %     end
    
    numSpks    = size(spks,1);
    iu         = allSpks(su).iu;
    bidsID     = allSpks(su).bidsID;
    sesh       = allSpks(su).sesh;
    wirename   = allSpks(su).wirename;
    bundlename = allSpks(su).bundlename;
    
    lfpName = [allMicro, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'];
    
    %% ONLY LOAD IN LFP IF IT IS A NEW ONE
    if ~strcmp(lfpName, prevLFPname)
        load(lfpName, 'data');
        prevLFPname = lfpName;
    end
    
    %% SELECT BUNDLE LFP
    selChan      = contains(data.label, bundlename);
    cfg          = [];
    cfg.channel  = data.label(selChan);       % all 8 channels
    bundLFP      = ft_selectdata(cfg, data);  % select
    
    numChan      = size(bundLFP.label,1);
    
%     %% WHICH TRIALS DOES THIS SU INDEX?
%     idxTrl = allSpks(su).idxTrl;
%     
%     %% IS THIS A (p)IU?
%     isIU = allSpks(su).iu>0;

    %% ENCODING AND RETRIEVAL TRIGGER
    encTrig  = round(allSpks(su).encTrigger(allSpks(su).hitsIdx,[1 3])*1000);
    retTrig  = round(allSpks(su).retTrigger(allSpks(su).hitsIdx,[1 3])*1000);
    
    %% EXTENDING TRIALS FOR FREQ-ANAL (to get the phase values without nan's at the sides)
    encTrigEx = [encTrig(:,1)-5000 encTrig(:,2)+5000];
    retTrigEx = [retTrig(:,1)-5000 retTrig(:,2)+5000];
    
    %% REDEFINE TRIALS ACCORDING TO encTrig (for AR)
    cfg        = [ ];
    cfg.trl    = [encTrigEx(:,1) encTrigEx(:,2) zeros(size(encTrigEx,1),1)];
    encLFPex   = ft_redefinetrial(cfg, bundLFP);
    
    cfg        = [ ];
    cfg.trl    = [encTrig(:,1) encTrig(:,2) zeros(size(encTrig,1),1)];
    encLFP     = ft_redefinetrial(cfg, bundLFP); % for AR
    
    %% DEMEAN AND ORTHOGONALIZE ENC
    cfg        = [];
    cfg.demean = 'yes';
    encLFPex   = ft_preprocessing(cfg, encLFPex); % demean
    encLFPex   = orthogonVec(encLFPex);           % orthogonalize
    
    %     encLFP     = ft_preprocessing(cfg, encLFP);
    
    %% REDEFINE TRIALS ACCORDING TO retTrig
    cfg        = [ ];
    cfg.trl    = [retTrigEx(:,1) retTrigEx(:,2) zeros(size(retTrigEx,1),1)];
    retLFPex   = ft_redefinetrial(cfg, bundLFP);
    
    cfg        = [ ];
    cfg.trl    = [retTrig(:,1) retTrig(:,2) zeros(size(retTrig,1),1)];
    retLFP     = ft_redefinetrial(cfg, bundLFP);
    
    %% DEMEAN AND ORTHOGONALIZE RET
    cfg        = [];
    cfg.demean = 'yes';
    retLFPex   = ft_preprocessing(cfg, retLFPex);
    retLFPex   = orthogonVec(retLFPex);
    
    %     retLFP     = ft_preprocessing(cfg, retLFP);
    
    %% LOOP OVER TRIALS
    for trl = 1 : size(encTrig,1)
        
        %% EXTRACT FREQSPEC
        cfgtf        = [];
        cfgtf.method = 'wavelet';
        cfgtf.width  = 6;
        cfgtf.toi    = 'all';
        %         cfgtf.foi    = logspace(log10(1),log10(140));
        cfgtf.foi = 1:1:45;
        cfgtf.output = 'fourier';
        cfgtf.trials = trl;
        
        %% LOW
        fspecEncLow        = ft_freqanalysis(cfgtf, encLFPex);
        fspecRetLow        = ft_freqanalysis(cfgtf, retLFPex);
        
        %% LOW-ENCODING
        encPhsLow = squeeze(fspecEncLow.fourierspctrm); % FOURIER SPECTRUM OF THAT TRIAL
        encPhsLow = encPhsLow(:,:,5001:end-5000); % SNIP WINGS
        encPhsLow = angle(encPhsLow); % TAKE ANGLE
        
        %% LOW-RETRIEVAL
        retPhsLow = squeeze(fspecRetLow.fourierspctrm); % FOURIER SPECTRUM OF THAT TRIAL
        retPhsLow = retPhsLow(:,:,5001:end-5000); % SNIP WINGS
        retPhsLow = angle(retPhsLow); % TAKE ANGLE
        
        %% HIGH
        cfgtf.width   = 12;
        cfgtf.foi     = 70:1:150;
        fspecEncHigh  = ft_freqanalysis(cfgtf, encLFPex);
        fspecRetHigh  = ft_freqanalysis(cfgtf, retLFPex);
        
        %% HIGH-ENCODING
        encPhsHigh = squeeze(fspecEncHigh.fourierspctrm); % FOURIER SPECTRUM OF THAT TRIAL
        encPhsHigh = encPhsHigh(:,:,5001:end-5000); % SNIP WINGS
        encPhsHigh = angle(encPhsHigh); % TAKE ANGLE
        
        %% HIGH-RETRIEVAL
        retPhsHigh = squeeze(fspecRetHigh.fourierspctrm); % FOURIER SPECTRUM OF THAT TRIAL
        retPhsHigh = retPhsHigh(:,:,5001:end-5000); % SNIP WINGS
        retPhsHigh = angle(retPhsHigh); % TAKE ANGLE
        
        %% SPIKES
        % WITHIN ENCODING TRIAL
        encSpkTrl = spks(spks>=encTrig(trl, 1) & spks <= encTrig(trl, 2))-encTrig(trl,1)+1; % ALL SPIKES THAT OCCUR WITHIN THE CURRENT TRIAL
        
        % WITHIN RETRIEVAL TRIAL
        retSpkTrl = spks(spks>=retTrig(trl, 1) & spks <= retTrig(trl, 2))-retTrig(trl,1)+1;
        
        %% PLUG SPKS INTO LFP PHASE
        %% LOW
        % ENCODING
        isAR                         = iqrAR(encLFP.trial{trl}, mkPlots);
        isAR                         = isAR(:,1:45,:);
        encPhsLow(isAR)              = NaN;
        phsDat                       = encPhsLow(:,:,encSpkTrl);
        allSpks(su).encPhsLow(trl,1) = {phsDat};           % in a cell because variable number of spikes per encoding trial
        
        % RETRIEVAL
        isAR                         = iqrAR(retLFP.trial{trl}, mkPlots);
        isAR                         = isAR(:,1:45,:);
        retPhsLow(isAR)              = NaN;
        phsDat                       = retPhsLow(:,:,retSpkTrl);
        allSpks(su).retPhsLow(trl,1) = {phsDat};           % in a cell because variable number of spikes per encoding trial
        
        %% HIGH
        % ENCODING
        isAR                          = iqrAR(encLFP.trial{trl}, mkPlots);
        isAR                          = isAR(:,1:81,:);
        encPhsHigh(isAR)              = NaN;
        phsDat                        = encPhsHigh(:,:,encSpkTrl);
        allSpks(su).encPhsHigh(trl,1) = {phsDat};           % in a cell because variable number of spikes per encoding trial
        
        % RETRIEVAL
        isAR                          = iqrAR(retLFP.trial{trl}, mkPlots);
        isAR                          = isAR(:,1:81,:);
        retPhsHigh(isAR)              = NaN;
        phsDat                        = retPhsHigh(:,:,retSpkTrl);
        allSpks(su).retPhsHigh(trl,1) = {phsDat};           % in a cell because variable number of spikes per encoding trial
        
    end
    
    %% WHICH CHANNEL IS THE CHOSEN ONE - LOW EDITION?
    % this is empty if there was not a single spike in any of the encoding/retrieval trials
    chanPhsEncLow = cat(3, allSpks(su).encPhsLow{:}); % every trial is {chan x freq x spk} and pooled over spk here
    chanPhsRetLow = cat(3, allSpks(su).retPhsLow{:});
    
    for chan = 1:8
        for freq = 1:size(chanPhsEncLow,2)
            % ENCODING
            extrPhs = squeeze(chanPhsEncLow(chan, freq,:));
            extrPhs(isnan(extrPhs)) = [];
            switch isempty(extrPhs)
                case 0
                    [~, zlvlEncLow(chan,freq)] = circ_rtest(extrPhs);
                case 1 % all spikes occur during artefacts
                    zlvlEncLow(chan,freq) = NaN;
            end
            
            % RETRIEVAL
            extrPhs = squeeze(chanPhsRetLow(chan, freq,:));
            extrPhs(isnan(extrPhs)) = [];
            switch isempty(extrPhs)
                case 0
                    [~, zlvlRetLow(chan,freq)] = circ_rtest(extrPhs);
                case 1 % all spikes occur during artefacts
                    zlvlRetLow(chan,freq) = NaN;
            end
        end
    end
    
    %% WHICH CHANNEL IS THE CHOSEN ONE - HIGH EDITION?
    chanPhsEncHigh = cat(3, allSpks(su).encPhsHigh{:}); % every trial is {chan x freq x spk} and pooled over spk here
    chanPhsRetHigh = cat(3, allSpks(su).retPhsHigh{:});
    
    for chan = 1:8
        for freq = 1:size(chanPhsEncHigh,2)
            % ENCODING
            extrPhs = squeeze(chanPhsEncHigh(chan, freq,:));
            extrPhs(isnan(extrPhs)) = [];
            switch isempty(extrPhs)
                case 0
                    [~, zlvlEncHigh(chan,freq)] = circ_rtest(extrPhs);
                case 1 % all spikes occur during artefacts
                    zlvlEncHigh(chan,freq) = NaN;
            end
            
            % RETRIEVAL
            extrPhs = squeeze(chanPhsRetHigh(chan, freq,:));
            extrPhs(isnan(extrPhs)) = [];
            switch isempty(extrPhs)
                case 0
                    [~, zlvlRetHigh(chan,freq)] = circ_rtest(extrPhs);
                case 1 % all spikes occur during artefacts
                    zlvlRetHigh(chan,freq) = NaN;
            end
        end
    end   
    
    
    %     for chan = 1:8
    %         [plvl(chan,1) zlvl(chan,1)] = circ_rtest(chanPhsEnc);
    %         [plvl(chan,2) zlvl(chan,2)] = circ_rtest(chanPhsRet);
    %     end
    
    %% OLD TO FIND THE HIGHEST CHANNEL OVER ALL FREQUENCIES
    %     encRetMax     = max(zlvlEnc, zlvlRet); % take maximum z-value of encoding vs. retrieval
    %     chanMax       = max(encRetMax, [], 2); % for each channel, take the maximum value over the frequency spectrum
    %     [~, chosChan] = max(chanMax);
    %     allSpks(su).favChan = chosChan;
    
    %% FIND THE CHANNEL WITH THE HIGHEST INPUT PER FREQUENCY (see cfgtf.foi)
        %     [~, chosChan]       = max(encRetMax, [], 1); % OLD for each channel, take the maximum value over the frequency spectrum
    encRetMaxLow        = max(zlvlEncLow, zlvlRetLow); % take maximum z-value of encoding vs. retrieval
    chanMax       = max(encRetMaxLow, [], 2); % for each channel, take the maximum value over the frequency spectrum
    [~, chosChanLow] = max(chanMax);
    allSpks(su).favChanLow = chosChanLow;
    
    encRetMaxHigh = max(zlvlEncHigh, zlvlRetHigh);
    chanMax       = max(encRetMaxHigh, [], 2); % for each channel, take the maximum value over the frequency spectrum
    [~, chosChanHigh] = max(chanMax);
    allSpks(su).favChanHigh = chosChanHigh;
end

save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_phs', 'allSpks', '-v7.3'); % with encoding and retrieval phases
allSpks = rmfield(allSpks, 'encPhsLow');
allSpks = rmfield(allSpks, 'retPhsLow');
allSpks = rmfield(allSpks, 'encPhsHigh');
allSpks = rmfield(allSpks, 'retPhsHigh');
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ', 'allSpks'); % light version with just the favourite channel
