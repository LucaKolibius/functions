%% IDEA OF THIS SCRIPT:
%  ZOOM INTO ALL SPIKES THAT OCCUR ON A SPECIFIC BUNDLE
%  FIRST FOR THE LOW FREQUENCIES AND THEN FOR THE HIGH FREQUENCIES GO
%  THROUGH THE SAME PROCEDURE:

%  IF MULTIPLE SPIKES ON THAT BUNDLE HAVE THE SAME INPUT WIRE THE INDEXED
%  TRIALS ARE ALL COLLAPSED TO THE FIRST SPIKE. ALL OTHER SPIKES ARE NAN'ED
%  THIS CREATES A POINTER TO WHICH TRIALS ON THAT WIRE THAT PROVIDES THE
%  INPUT LFP WERE INDEXED AND WHICH WERENT. AT THE SAME TIME IT PREVENTS
%  LFPs TO BE COUNTED DOUBLE.

clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat')

allSpks = rmfield(allSpks, 'idxTrlSing');
allSpks = rmfield(allSpks, 'favChan');
allSpks = rmfield(allSpks, 'idxTrlSingHi');
allSpks = rmfield(allSpks, 'idxTrlSingLw');

% if ~isfield(allSpks, 'idxTrlSing')
%     allSpks(1).idxTrlSing = [];
% end

% for su = 1 : size(allSpks,2)
%     allSpks(su).favChanLow = round(unifrnd(1,8));
%     allSpks(su).favChanHigh = round(unifrnd(1,8));
% end

allBids = {allSpks.bidsID};
allBund = {allSpks.bundlename};
allSesh = {allSpks.sesh};

for su = 1 : size(allSpks,2)
    disp(su)
%     if isnan(allSpks(su).idxTrlSing)
%         continue
%     end
    bundlename = allSpks(su).bundlename;
    bidsID     = allSpks(su).bidsID;
    sesh       = allSpks(su).sesh;
    
    spkBlock   = strcmp(allBids, bidsID) & strcmp(allBund, bundlename) & strcmp(allSesh, sesh); % all the spikes on that bundle
    spkBlock   = find(spkBlock == 1);
    
    favChanLow    = [allSpks(spkBlock).favChanLow]; % favourite channels of all the single units in that bundle (low frequencies)
    favChanHigh   = [allSpks(spkBlock).favChanHigh];
    
    for chan = 1:8
        %% LOW
        idxChan   = find(favChanLow==chan); % go through all 8 possible favourite channels.
        
        if isempty(idxChan) % if that wire does not supply the input to any SU, continue with the next one
            continue
        end
        
        idxTrls   = vertcat(allSpks(spkBlock(idxChan)).idxTrl); % combined index score for this input LFP over all SU
        newIdxTrl = any(idxTrls,1);
        
        allSpks(spkBlock(idxChan(1))).idxTrlSingLw     = newIdxTrl;
        
        if size(idxChan,2) > 1
            for extr = 2:size(idxChan,2)
                allSpks(spkBlock(idxChan(extr))).idxTrlSingLw = NaN;
            end
        end
    end
    
    for chan = 1:8
        %% HIGH
        idxChan   = find(favChanHigh==chan); % go through all 8 possible favourite channels.
        
        
        if isempty(idxChan) % if that wire does not supply the input to any SU, continue with the next one
            continue
        end
        
        idxTrls   = vertcat(allSpks(spkBlock(idxChan)).idxTrl); % combined index score for this input LFP over all SU
        newIdxTrl = any(idxTrls,1);
        
        allSpks(spkBlock(idxChan(1))).idxTrlSingHi     = newIdxTrl;
        
        if size(idxChan,2) > 1 % if multiple spikes get their LFP input from the same channel, nan all others after the first
            for extr = 2:size(idxChan,2)
                allSpks(spkBlock(idxChan(extr))).idxTrlSingHi = NaN;
            end
        end
    end
end
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')

