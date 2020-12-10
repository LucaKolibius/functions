clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat')
allSpks(1).idxTrlSing = [];

% for su = 1 : size(allSpks,2)
%     allSpks(su).favChan = round(unifrnd(1,8));
% end

allBids = {allSpks.bidsID};
allBund = {allSpks.bundlename};
allSesh = {allSpks.sesh};

for su = 1 : size(allSpks,2)
    disp(su)
    if isnan(allSpks(su).idxTrlSing)
        continue
    end
    bundlename = allSpks(su).bundlename;
    bidsID     = allSpks(su).bidsID;
    sesh       = allSpks(su).sesh;
    
    spkBlock   = strcmp(allBids, bidsID) & strcmp(allBund, bundlename) & strcmp(allSesh, sesh); % all the spikes on that bundle
    spkBlock   = find(spkBlock == 1);
    
    favChan    = [allSpks(spkBlock).favChan]; % favourite channels of all the single units in that bundle
    
    for chan = 1:8
        idxChan   = find(favChan==chan); % go through all 8 possible favourite channels.
        
        if isempty(idxChan) % if that wire does not supply the input to any SU, continue with the next one
            continue
        end
        
        idxTrls   = vertcat(allSpks(spkBlock(idxChan)).idxTrl); % combined index score for this input LFP over all SU
        newIdxTrl = any(idxTrls,1);
        
        allSpks(spkBlock(idxChan(1))).idxTrlSing     = newIdxTrl;
        
        if size(idxChan,2) > 1
            for extr = 2:size(idxChan,2)
                allSpks(spkBlock(idxChan(extr))).idxTrlSing = NaN;
            end
        end
        
    end
end
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')

