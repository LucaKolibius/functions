% clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp.mat')
isCN = false(length(allSpks),1);
for su = 1: length(allSpks)
    lineLength = fprintf('%d SU out of %d SU done (%.2f%%).\n', su, length(allSpks), su/length(allSpks)*100);
    
    %% ignore SU for now
%     if allSpks(su).iu == 0
%         continue
%     end
    
    %% TRIGGER IS IN SECONDS
    
    % RETRIEVAL TRIGGER (CUE)
    trigRet = allSpks(su).retTrigger(allSpks(su).hitsIdx, [1]);
    % ENCODING TRIGGER (CUE)
    trigEnc = allSpks(su).encTrigger(allSpks(su).hitsIdx, [1]);
    % REINSTATED EPISODES
    idxTrl     = allSpks(su).idxTrl;
    
    % SPIKETIMES (IN SECONDS)
    clusterSpikes = allSpks(su).spks/1000;
    
    
    spksEncBL = insertSpikeNumber(trigEnc, clusterSpikes, [-1.000 -0.300]);
    spksRetBL = insertSpikeNumber(trigRet, clusterSpikes, [-1.000 -0.300]);
    spksBL    = vertcat(spksEncBL, spksRetBL);
    
    spksEncTr = insertSpikeNumber(trigEnc, clusterSpikes, [ 0.300  0.700]);
    spksRetTr = insertSpikeNumber(trigRet, clusterSpikes, [ 0.300  0.700]);
%     spksTr    = [spksEncTr spksRetTr]; % ENCODING AND RETRIEVAL POSTCUE ACTIVITY
    spksTr    = spksEncTr; %% ONLY ENCODING
    
    for trl = 1:size(trigEnc,1)
        
%     spksTr    = spksTr(idxTrl,:); % FOR ESN ANALYSIS
    spkTr = spksTr(trl,:);
    medFir    = median(spksTr,'all');
    
    blMean    = mean(spksBL);
    blStd     = std(spksBL);
    th        = blMean + blStd * 5;
    
    if medFir >= 2 && medFir > th
        isCN(su) = true;
    end
    
    end
    fprintf(repmat('\b', 1, lineLength));
end

sum(isCN)

isESN = [allSpks.iu] == 1 | [allSpks.iu] == 2;

sum(isESN(isCN))


%%
noCN = ~isCN
noCN = find(noCN)