%% SFC enc/ret X idx/ndx
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

allPhsDat.encNdx = [];
allPhsDat.encIdx = [];
allPhsDat.retNdx = [];
allPhsDat.retIdx = [];
for su = 1: 4%size(allSpks,2)
    if allSpks(su).iu == 0
        continue
    end
    
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
    favChan = allSpks(su).favChanHigh;

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
        cfgtf.toi    = 'all';
        cfgtf.channel = favChan;
        %         cfgtf.foi    = logspace(log10(1),log10(140));
        cfgtf.output = 'fourier';
        cfgtf.trials = trl;
        
        
        %% HIGH
        cfgtf.width   = 12;
        cfgtf.foi     = 60:2:160;
        fspecEncHigh  = ft_freqanalysis(cfgtf, encLFPex);
        fspecRetHigh  = ft_freqanalysis(cfgtf, retLFPex);
        
        %% HIGH-ENCODING
        encPhsHigh = squeeze(fspecEncHigh.fourierspctrm); % FOURIER SPECTRUM OF THAT TRIAL
        encPhsHigh = encPhsHigh(:,5001:end-5000); % SNIP WINGS
        encPhsHigh = angle(encPhsHigh); % TAKE ANGLE
        
        %% HIGH-RETRIEVAL
        retPhsHigh = squeeze(fspecRetHigh.fourierspctrm); % FOURIER SPECTRUM OF THAT TRIAL
        retPhsHigh = retPhsHigh(:,5001:end-5000); % SNIP WINGS
        retPhsHigh = angle(retPhsHigh); % TAKE ANGLE
        
        %% SPIKES
        % WITHIN ENCODING TRIAL
        encSpkTrl = spks(spks>=encTrig(trl, 1) & spks <= encTrig(trl, 2))-encTrig(trl,1)+1; % ALL SPIKES THAT OCCUR WITHIN THE CURRENT TRIAL
        
        % WITHIN RETRIEVAL TRIAL
        retSpkTrl = spks(spks>=retTrig(trl, 1) & spks <= retTrig(trl, 2))-retTrig(trl,1)+1;
        
        %% PLUG SPKS INTO LFP PHASE        
        %% HIGH
        % ENCODING
        isAR                          = iqrAR(encLFP.trial{trl}, mkPlots);
        isAR                          = squeeze(isAR(favChan,1:length(cfgtf.foi),:));
        encPhsHigh(isAR)              = NaN;
        phsDat                        = encPhsHigh(:,encSpkTrl);
        
        switch idxTrl(trl)
            case 0
                allPhsDat(su).encNdx = [allPhsDat.encNdx, phsDat];
            case 1
                allPhsDat(su).encIdx = [allPhsDat.encIdx, phsDat];
        end
                   
        % RETRIEVAL
        isAR                          = iqrAR(retLFP.trial{trl}, mkPlots);
        isAR                          = squeeze(isAR(favChan,length(cfgtf.foi),:));
        retPhsHigh(isAR)              = NaN;
        phsDat                        = retPhsHigh(:,retSpkTrl);

        switch idxTrl(trl)
            case 0
                allPhsDat(su).retNdx = [allPhsDat.retNdx, phsDat];
            case 1
                allPhsDat(su).retIdx = [allPhsDat.retIdx, phsDat];
        end        
    end
    
   
end

save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allPhsDat.mat', 'allPhsDat');

% figure(1); clf;
% subplot(211); hold on;
% hand = histogram(allPhsDat.encIdx,'Normalization', 'probability');
% hand.FaceColor = [0.8510    0.3255    0.0980];
% hand = histogram(allPhsDat.encNdx,'Normalization', 'probability');
% hand.FaceColor = [0.0980 0.8510    0.3255];
% 
% subplot(212); hold on;
% hand = histogram(allPhsDat.retIdx,'Normalization', 'probability');
% hand.FaceColor = [0.8510    0.3255    0.0980];
% hand = histogram(allPhsDat.retNdx,'Normalization', 'probability');
% hand.FaceColor = [0.0980 0.8510    0.3255];
