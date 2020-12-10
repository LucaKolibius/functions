function preCuePowDiff
lfpDir = dir('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_noSPKINT.mat');  % NO SPK INT
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')

skippedDat = [];
mtrack = [];

allSUPow.idx  = [];
allSUPow.ndx  = [];

for spk = 1 : length(allSpks)
    
    if isnan(allSpks(spk).idxTrlSing)
        continue
    end
    
    if isempty(allSpks(spk).favChan)
        error('spPhase_encRet fully completed!')
    end
    
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
    curBund = allSpks(spk).bundlename;
    
    doneLoad = 0;
    while doneLoad == 0
        try
            curWire = [curBund, num2str(allSpks(spk).favChan)];
            doneLoad = 1;
        catch
            disp('Trying to extract favourite channel again in 60s');
            pause(60);
        end
    end
    
    idxTrl  = allSpks(spk).idxTrlSing;     % WHICH TRIALS DOES THAT SU INDEX?
    encTrig  = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1])*1000);

    %% LOAD IN THE LFP-DATA
    try
        load([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    catch
        skippedDat = [skippedDat; {bidsID} {sesh}];
        disp('skipDat');
        error('No files should be skipped anymore.')
    end
    
   
    %% REDEFINE TRIALS ACCORDING TO encTrig
    cfg          = [ ];
    cfg.trl      = [encTrig-1000 encTrig-1 zeros(size(encTrig))];
    microLFP     = ft_redefinetrial(cfg, data);
      
    
    %% ONLY DEMEAN & ORTHOGONALIZE
        % normalize LFP variance
    %     data.trial = {(data.trial{1} - mean(data.trial{1},2)) ./ std(data.trial{1},0,2)};

    
    cfg        = [];
    cfg.demean = 'yes';
    microLFP   = ft_preprocessing(cfg, microLFP);
    microLFP   = orthogonVec(microLFP);           % orthogonalize
    
    
    %% SELECT CHANNEL THAT REFLECTS SU INPUT (see spkPhase_encRet)
    cfg          = [];
    cfg.channel  = curWire;
    microLFP     = ft_selectdata(cfg, microLFP); % select
    
    idxTrlPow = [];
    ndxTrlPow = [];
    goodTrl   = ones(1, length(idxTrl));
    for trl = 1 : length(encTrig)
      
      %% AR
      isAR = iqrAR(microLFP.trial{trl},0);
      isAR = squeeze(isAR(1,1,:));

      if any(isAR)
          goodTrl(trl) = 0;
          continue
      end
 
        % FIELDTRIP POWER
        cfg = [];
        cfg.output    = 'pow';
        cfg.channel   = 'all';
        cfg.method    = 'mtmfft';
        cfg.pad       = 4; % padding to increase frequency resolution
        cfg.foi       = [1:0.25:200];
        cfg.taper     = 'hanning';
        cfg.trials    = trl;
        trlPow        = ft_freqanalysis(cfg, microLFP);
        
        switch idxTrl(trl)
            case 0
                ndxTrlPow = [ndxTrlPow; trlPow.powspctrm];
            case 1
                idxTrlPow = [idxTrlPow; trlPow.powspctrm];
                
        end
        
    end
    
    if sum(idxTrl(logical(goodTrl))) > 0
        
        allSUPow.idx  = [ allSUPow.idx  {idxTrlPow}  ];
        allSUPow.ndx  = [ allSUPow.ndx  {ndxTrlPow}  ];
        
        addthis = [spk; length(allSUPow.idx)];
        mtrack  = [mtrack, addthis];
    end
    
    %% UPDATE PROGRESS
    lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).', spk, length(allSpks), spk/length(allSpks)*100);
%     fprintf(repmat('\b',1,lineLength))

end
hz = trlPow.freq;
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\preCuePowDiff_orthDeMea.mat', 'allSUPow', 'hz')

end % END OF FUNCTION