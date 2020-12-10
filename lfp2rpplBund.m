%% NEW VERSION FOR MICRO LFP WITHOUT SPKINT THAT ONLY CONSIDERES BUNDLES IN WHICH I HAVE HIPPOCAMPAL UNITS
clear
addpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\functions');
addpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\toolboxes\fieldtrip-20200603')

lfpDir = dir('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_noSPKINT.mat');  % CHANGED TO NO SPK INT
% artFold = 'Z:\hanslmas-ieeg-compute\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS
% artFold = 'X:\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')

allSUPow.idx  = [];
allSUPow.ndx  = [];
allSUPow.dff  = [];
allFreqRes    = [];

for spk = 1 : length(allSpks)
    
    if any(isnan(allSpks(spk).idxTrlSing))
        continue
    end
    
    if isempty(allSpks(spk).favChan)
        continue
    end
    
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
%     allBund = {allSpks.bundlename};
    curBund = allSpks(spk).bundlename;
%     curWire = [curBund, num2str(allSpks(spk).favChan)];
    favChan = allSpks(spk).favChan(45:50);
    idxTrl  = allSpks(spk).idxTrlSing;     % WHICH TRIALS DOES THAT SU INDEX?
    encTrig = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1 3])*1000);

    %% LOAD IN THE LFP-DATA
    % There is no indendent LFP for sub-1007_S1b (it is in the same
    % recording as S1, it's just another session)
    load([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    
    %% REDEFINE TRIALS ACCORDING TO encTrig
    cfg          = [ ];
    cfg.trl      = [encTrig(:,1) encTrig(:,2) zeros(size(encTrig,1))];
    microLFP     = ft_redefinetrial(cfg, data);
    microLFP     = rmfield(microLFP, 'trialinfo');
    
    %% DEMEAN & ORTHOGONALIZE   
    cfg        = [];
    cfg.demean = 'yes';
    microLFP   = ft_preprocessing(cfg, microLFP);
    microLFP   = orthogonVec(microLFP);           % orthogonalize
    
    %% SELECT CHANNEL THAT REFLECTS SU INPUT (see spkPhase_encRet)
%     cfg          = [];
%     cfg.channel  = curWire;
%     microLFP     = ft_selectdata(cfg, microLFP); % select

%% THIS CHANNEL SELECTION COULD BE A DIFFERENT WIRE FOR EACH FREQUENCY
allBund = cellfun(@(x) x(1:end-1), microLFP.label, 'un', 0);
bundIdx = strcmp(allBund, curBund);

cfg         = [];
cfg.channel = microLFP.label(bundIdx);
microLFP    = ft_selectdata(cfg, microLFP);
%     %% DETECT RIPPLES
%     [rpplsWire, rpplLen, bndLab] = calcRppl (microLFP);
%     allSpks(spk).rppls   = rpplsWire;
%     allSpks(spk).rpplLen = rpplLen;
    
    %% CALCULATE TRIAL RIPPLE POWER (80-140hz)
    [idxTrlPow, ndxTrlPow, goodTrl, freqRes] = calcRipplPow (microLFP, idxTrl, encTrig, favChan);
    
if sum(idxTrl(logical(goodTrl))) > 0
    
    allSUPow.idx  = [ allSUPow.idx  {idxTrlPow}  ];
    allSUPow.ndx  = [ allSUPow.ndx  {ndxTrlPow}  ];
    allSUPow.dff  = [ allSUPow.dff {}];
    
end

allFreqRes = [allFreqRes; freqRes];
end
% save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_rppls.mat', 'allSpks');

save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\rpplPowDiff_orthDeMea.mat', 'allSUPow', 'allFreqRes')
