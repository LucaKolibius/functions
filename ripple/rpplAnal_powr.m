function rpplAnal_powr

whereAmI(0)
global prePath;
addpath([prePath, 'Luca\functions']);
addpath([prePath, 'Luca\toolboxes\fieldtrip-20200310']); ft_defaults;
addpath([prePath, 'Luca\toolboxes\fieldtrip-20200603']); ft_defaults;
rmpath(genpath('\\analyse4.psy.gla.ac.uk\project0310\RDS\Common\mcode'));
% lfpDir = dir([prePath, 'Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_noSPKINT.mat']);  % CHANGED TO NO SPK INT
lfpDir = dir([prePath, 'Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_*.mat']);  % spkInt in laptop and noSpkInt on local
load([prePath, 'Luca\data\allSbj\allSpksHZ.mat'], 'allSpks')

allSUPow.idx  = [];
allSUPow.ndx  = [];
allSUPow.dff  = [];
allFreqRes    = [];

for spk = 1 : length(allSpks)
    
        % PREVENT DOUBLES
        if any(isnan(allSpks(spk).idxTrlSingHi))
            continue
        end
    
        % NO IDX TRIALS
        if sum(allSpks(spk).idxTrlSingHi) == 0
            continue
        end
        
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
    %     allBund = {allSpks.bundlename};
    curBund = allSpks(spk).bundlename;
    %     curWire = [curBund, num2str(allSpks(spk).favChan)];
    %     favChan = allSpks(spk).favChan(45:50);
    idxTrl  = allSpks(spk).idxTrlSingHi;     % WHICH TRIALS DOES THAT WIRE INDEX?
    encTrig = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1 3])*1000);
    favChan = allSpks(spk).favChanHigh;
    
    %% LOAD IN THE LFP-DATA
    % There is no indendent LFP for sub-1007_S1b (it is in the same
    % recording as S1, it's just another session)
    load([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    
    % LAPTOP AND NO SPK INT
    %         abc = dir([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_*.mat']);
    %         load([abc.folder, filesep, abc.name], 'data');
    
    %% REDEFINE TRIALS ACCORDING TO encTrig
    cfg          = [ ];
    cfg.trl      = [encTrig(:,1)-100 encTrig(:,2)+100 zeros(size(encTrig,1))];
    microLFPex = ft_redefinetrial(cfg, data);
    microLFPex = rmfield(microLFPex, 'trialinfo');
    
    cfg          = [ ];
    cfg.trl      = [encTrig(:,1) encTrig(:,2) zeros(size(encTrig,1))];
    microLFP     = ft_redefinetrial(cfg, data);
    microLFP     = rmfield(microLFP, 'trialinfo');
    
    %% DEMEAN & ORTHOGONALIZE
    cfg        = [];
    cfg.demean = 'yes';
    microLFP   = ft_preprocessing(cfg, microLFP);
    microLFP   = orthogonVec(microLFP);           % orthogonalize
    
    cfg        = [];
    cfg.demean = 'yes';
    microLFPex   = ft_preprocessing(cfg, microLFPex);
    microLFPex   = orthogonVec(microLFPex);           % orthogonalize

    %% SELECT CHANNEL THAT REFLECTS SU INPUT (see spkPhase_encRet)
%         cfg          = [];
%         cfg.channel  = curWire;
%         microLFP     = ft_selectdata(cfg, microLFP); % select
    
    %% CHANNEL SELECTION COULD BE A DIFFERENT WIRE FOR EACH FREQUENCY
    allBund = cellfun(@(x) x(1:end-1), microLFP.label, 'un', 0);
    bundIdx = strcmp(allBund, curBund);
    
    cfg         = [];
    cfg.channel = microLFP.label(bundIdx);
    microLFP    = ft_selectdata(cfg, microLFP);
    microLFPex    = ft_selectdata(cfg, microLFPex);

    %% DETECT RIPPLES
        %% CALCULATE TRIAL RIPPLE POWER (70-150hz)
        [idxTrlPow, ndxTrlPow, freqRes] = calcRipplPow (microLFP, microLFPex, idxTrl, encTrig, favChan);
    
    
            allSUPow.idx  = [ allSUPow.idx  {idxTrlPow}  ];
            allSUPow.ndx  = [ allSUPow.ndx  {ndxTrlPow}  ];
    
    
            allFreqRes = [allFreqRes; freqRes];
            
            %% UPDATE PROGRESS
            lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
            lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
            lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
            lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
            lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', spk, length(allSpks), spk/length(allSpks)*100);
            
            
end
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_rppls_pow.mat', 'allSUPow', 'allFreqRes', '-v7.3');
end % OF FUNCTION