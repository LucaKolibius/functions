data;

%% 2ms before 6ms after
clear; clc;
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHz.mat', 'allSpks')
cd('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\microLFP\');
allRec = dir('sub-*noSPKINT.mat');

allBids = {allSpks.bidsID};
allSesh = {allSpks.sesh};
for rec = 1:length(allRec)
    load([allRec(1).folder, filesep, allRec(rec).name], 'data');
    curTxt = allRec(rec).name;
    bids = curTxt(1:8);
    sesh = curTxt(10:11);
    
    spkIdx = strcmp(bids, allBids) & strcmp(sesh, allSesh);
    spkIdx = find(spkIdx);
    numSpks = length(spkIdx);
    % loop through spikes
    for su = 1 : numSpks
        curWire = allSpks(spkIdx(su)).wirename;
        wireIdx = strcmp(data.label, curWire);
        spks = round(allSpks(spkIdx(su)).spks);
        
        for spk = 1:length(spks) % Loop through all individual spikes
            spkTime = spks(spk);
            intWin = [2 6];
            while spkTime-intWin(1) < 1
                intWin(1) = intWin(1)-1;
            end
            
            while spkTime+intWin(2) > size(data.trial{1},2)
                intWin(2) = intWin(2)-1;
            end
            
%             hold off
%             plot(data.trial{1}(wireIdx, spkTime-20:spkTime+20)); hold on;
            startVal = data.trial{1}(wireIdx,spkTime- intWin(1));
            endVal   = data.trial{1}(wireIdx,spkTime+ intWin(2));
            
            intVals  = linspace(startVal, endVal, sum(intWin)+1);
            
            data.trial{1}(wireIdx, spkTime -intWin(1)+1:spkTime+ (intWin(2)-1)) = intVals(2:end-1);
            
%             plot(data.trial{1}(wireIdx,spkTime-20:spkTime+20))
%             pause(1);
        end
    end
    
    save(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data\microLFP\mySpkInt\',bids, '_', sesh, '_onlyMicroLFP_RAW_1000DS_mSPKint.mat'], 'data');
end