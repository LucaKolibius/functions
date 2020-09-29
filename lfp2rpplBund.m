%% NEW VERSION FOR MICRO LFP WITHOUT SPKINT THAT ONLY CONSIDERES BUNDLES IN WHICH I HAVE HIPPOCAMPAL UNITS
clear
lfpDir = dir('X:\Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_noSPKINT.mat');  % CHANGED TO NO SPK INT
% artFold = 'Z:\hanslmas-ieeg-compute\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS
artFold = 'X:\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS
load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')
counter = 1;
skippedDat = [];

for spk = 1 : length(allSpks)
    
    %% SKIP REPEATING MICROWIRE BUNDLES
    if spk > 1
        if and(and(contains(allSpks(spk-1).bidsID, allSpks(spk).bidsID), strcmp(allSpks(spk-1).sesh, allSpks(spk).sesh) ), strcmp(allSpks(spk).bundlename, allSpks(spk-1).bundlename)); % same subject + session
            continue
        end
    end
    
    %% GET: bidsID + sesh
    bidsID = allSpks(spk).bidsID;
    sesh   = allSpks(spk).sesh;
    allBund = {allSpks.bundlename};
    curBund = allSpks(spk).bundlename;
    
    %% LOAD IN THE LFP-DATA
    % There is no indendent LFP for sub-1007_S1b (it is in the same
    % recording as S1, it's just another session)
    try
        load([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    catch
        error('You should not have to skip any more data (all is fixed)')
        %         skippedDat = [skippedDat; {bidsID} {sesh}];
        %         disp('skipDat');
        %         continue
    end
    
    
    %% SELECT BUNDLE OF INTEREST
    dataBundLab = cellfun(@(x) x(1:end-1) , data.label, 'un', 0);
    loadBunds   = strcmp(curBund, dataBundLab);
    
    cfg          = [];
    cfg.channel  = data.label(loadBunds);       % all 8 channels
    data         = ft_selectdata(cfg, data);  % select
    
    %% DETECT RIPPLES
    [~, rpplsWire, bndLab] = calcRppl (data, bidsID, sesh, artFold, 1);
    
    for bund = 1 : size(rpplsWire,1)
        
        %% WHICH TRIALS IN THAT BUNDLE ARE INDEXED?
        sameBund = and(and(contains({allSpks.bidsID}, bidsID), strcmp({allSpks.sesh}, sesh) ), contains(allBund, bndLab(bund))); % same subject + session + bundle
        idxTrl   = any(vertcat(allSpks(sameBund).idxTrl),1); % these are the trials that are indexed in that wire
        
        %% SAVE INFO TO NEW VARIABLE
        rpplBund(counter).bidsID     = bidsID;
        rpplBund(counter).sesh       = sesh;
        rpplBund(counter).bundname   = bndLab{bund};
        rpplBund(counter).rppls      = rpplsWire{bund};
        rpplBund(counter).idxTrlBund = idxTrl;
        rpplBund(counter).encTrigger = allSpks(spk).encTrigger(allSpks(spk).hitsIdx,:);
        
        counter  = counter + 1;
 
    end
end
save('X:\Luca\data\allSbj\rpplBund_noSpkIntOrth.mat', 'rpplBund', 'skippedDat');