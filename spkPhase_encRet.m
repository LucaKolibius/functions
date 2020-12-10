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
for su = 1 : size(allSpks,2)
    
    disp(su)
    
    if ~isempty(allSpks(su).favChan)
        continue
    end
    
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
    
    %% WHICH TRIALS DOES THIS SU INDEX?
    idxTrl = allSpks(su).idxTrl;
    
    %% IS THIS A (p)IU?
    isIU = allSpks(su).iu>0;
    
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
        cfgtf.width  = 8;
        cfgtf.toi    = 'all';
        cfgtf.foi    = logspace(log10(1),log10(140));
        cfgtf.output = 'fourier';
        cfgtf.trials = trl;
        
        fspecEnc        = ft_freqanalysis(cfgtf, encLFPex);
        fspecRet        = ft_freqanalysis(cfgtf, retLFPex);
        
        % ENCODING
        encPhs = squeeze(fspecEnc.fourierspctrm); % FOURIER SPECTRUM OF THAT TRIAL
        encPhs = encPhs(:,:,5001:end-5000); % SNIP WINGS
        encPhs = angle(encPhs); % TAKE ANGLE
        
        % RETRIEVAL
        retPhs = squeeze(fspecRet.fourierspctrm); % FOURIER SPECTRUM OF THAT TRIAL
        retPhs = retPhs(:,:,5001:end-5000); % SNIP WINGS
        retPhs = angle(retPhs); % TAKE ANGLE
        
        %% SPIKES
        % WITHIN ENCODING TRIAL
        encSpkTrl = spks(spks>=encTrig(trl, 1) & spks <= encTrig(trl, 2))-encTrig(trl,1)+1; % ALL SPIKES THAT OCCUR WITHIN THE CURRENT TRIAL
        
        % WITHIN RETRIEVAL TRIAL
        retSpkTrl = spks(spks>=retTrig(trl, 1) & spks <= retTrig(trl, 2))-retTrig(trl,1)+1;
        
        %% PLUG SPKS INTO LFP PHASE
        % ENCODING
        isAR                      = iqrAR(encLFP.trial{trl}, mkPlots);
        encPhs(isAR)              = NaN;
        phsDat                    = encPhs(:,:,encSpkTrl);
        allSpks(su).encPhs(trl,1) = {phsDat};           % in a cell because variable number of spikes per encoding trial
        
        % RETRIEVAL
        isAR                      = iqrAR(retLFP.trial{trl}, mkPlots);
        retPhs(isAR)              = NaN;
        phsDat                    = retPhs(:,:,retSpkTrl);
        allSpks(su).retPhs(trl,1) = {phsDat};           % in a cell because variable number of spikes per encoding trial
        
        
    end
    
    %% WHICH CHANNEL IS THE CHOSEN ONE?
    % this is empty if there was not a single spike in any of the encoding/retrieval trials
    chanPhsEnc = cat(3, allSpks(su).encPhs{:}); % every trial is {chan x freq x spk} and pooled over spk here
    chanPhsRet = cat(3, allSpks(su).retPhs{:});
    
    for chan = 1:8
        for freq = 1:size(chanPhsEnc,2)
            % ENCODING
            extrPhs = squeeze(chanPhsEnc(chan, freq,:));
            extrPhs(isnan(extrPhs)) = [];
            switch isempty(extrPhs)
                case 0
                    [~, zlvlEnc(chan,freq)] = circ_rtest(extrPhs);
                case 1 % all spikes occur during artefacts
                    zlvlEnc(chan,freq) = NaN;
            end
            
            % RETRIEVAL
            extrPhs = squeeze(chanPhsRet(chan, freq,:));
            extrPhs(isnan(extrPhs)) = [];
            switch isempty(extrPhs)
                case 0
                    [~, zlvlRet(chan,freq)] = circ_rtest(extrPhs);
                case 1 % all spikes occur during artefacts
                    zlvlRet(chan,freq) = NaN;
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
    encRetMax           = max(zlvlEnc, zlvlRet); % take maximum z-value of encoding vs. retrieval
    [~, chosChan]       = max(encRetMax, [], 1); % for each channel, take the maximum value over the frequency spectrum
    allSpks(su).favChan = chosChan;
    
end

save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_phs', 'allSpks', '-v7.3'); % with encoding and retrieval phases
allSpks = rmfield(allSpks, 'encPhs');
allSpks = rmfield(allSpks, 'retPhs');
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ', 'allSpks'); % light version with just the favourite channel
