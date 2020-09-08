clear all

load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')
counter = 1;

for spk = 1 : length(allSpks)
    
    if spk > 1
        if strcmp(allSpks(spk).bidsID, allSpks(spk-1).bidsID) & strcmp(allSpks(spk).sesh, allSpks(spk-1).sesh) % same sesh and same sbj
            continue
        end
    end
    
    % SAVE INFO ABOUT SUB & SESH
    supps{counter, 1} = allSpks(spk).bidsID;
    supps{counter, 2} = allSpks(spk).sesh;
    
    % CD TO ELEC LOC
    subj = sub_ID_conversion(allSpks(spk).bidsID, 1);
    
    if strcmp(subj, 'P7_ERL')
        subj = 'P07ERL';
    end
    
    cd(['Z:\hanslmas-ieeg-compute\Luca\data\', subj])
    cd(allSpks(spk).sesh)
    abc = dir('20*'); cd(abc(1).name);
    cd('advancedAnalysis\')
    cd('elecLoc')
    
    % LOAD ELEC LOC
    abc = dir('elecLoc*');
    load(abc.name)
    
    % SAVE INFO ABOUT BUNDLES
    supps{counter, 3} = size(elecLoc,1); % NUMBER OF BUNDLES
    supps{counter, 4} = sum(strcmp(elecLoc(:,2), 'hipp')); % NUMBER OF BUNDLES IN HIPP
    
    % SAVE INFO ABOUT TRIAL NUMBER AND HITS
    supps{counter, 5} = size(allSpks(spk).encTrigger, 1);
    supps{counter, 6} = size(allSpks(spk).hitsIdx, 1);
    
    % NUMBER OF BUNDLES WITH SU
    sameSesh = and(contains({allSpks.bidsID}, allSpks(spk).bidsID), strcmp({allSpks.sesh}, allSpks(spk).sesh) ); % same subject + session
    
    allBund = {allSpks(sameSesh).bundlename};
    allBund = unique(allBund);
    
    supps{counter, 7} = length(allBund);
    
    % NUMBER OF SU
    supps{counter,8} = sum(sameSesh);
    
    % NUMBER OF PUTATIVE INDEX UNITS AND INDEX UNITS
    supps{counter, 9}  = sum([allSpks(sameSesh).iu] == 2);
    supps{counter, 10} = sum([allSpks(sameSesh).iu] == 1);
    
    counter = counter + 1;
end

sbj = unique([supps(:,1)]);

mean([supps{:,6}]./ [supps{:,5}])
std([supps{:,6}]./ [supps{:,5}])
