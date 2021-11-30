clear 
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp_noCN.mat', 'allSpks');counter = 1;

skipit = zeros(size(allSpks,2),1);
for spk = 1 : length(allSpks)
    
    if skipit(spk) == 1
        continue
    end

    
%     if spk > 1
%         if strcmp(allSpks(spk).bidsID, allSpks(spk-1).bidsID) & strcmp(allSpks(spk).sesh, allSpks(spk-1).sesh) % same sesh and same sbj
%             continue
%         end
%     end
    
    % SAVE INFO ABOUT SUB & SESH
    supps{counter, 1} = allSpks(spk).bidsID;
    supps{counter, 2} = allSpks(spk).sesh;
    
    % CD TO ELEC LOC
    subj = sub_ID_conversion(allSpks(spk).bidsID, 1);
    
    if strcmp(subj, 'P7_ERL')
        subj = 'P07ERL';
    end
    
    cd(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data\', subj])
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
    
    % HOW MANY TRIALS ARE REINSTATED FOR ESNs?
    ESNidx = ([allSpks.iu] == 1 | [allSpks.iu] == 2) & sameSesh == 1;
    numReinst = vertcat(allSpks(ESNidx).idxTrl);
    numReinst = sum(numReinst,2);
    stdReinst = std(numReinst); %??
    numReinst = mean(numReinst,1);
    supps{counter, 11} = vertcat(allSpks(sameSesh).idxTrl);
    
    counter = counter + 1;
        
    % SKIP ALL OTHER SU IN THAT SUBJECT AND SESSION
    skipit(sameSesh) = 1;
end

sbj = unique([supps(:,1)]);

mean([supps{:,6}]./ [supps{:,5}])
std([supps{:,6}]./ [supps{:,5}])


%% hitrate
hrate = ([supps{:,6}]./ [supps{:,5}]);

allSub = supps(:,1);
allSes = supps(:,2);
allHit = [];

while ~isempty(allSub)
    idx = strcmp(allSub, allSub(1)); % all sbj & strcmp(allSes, supps(sesh,2))
    
    allHit = [allHit; mean(hrate(idx))];
    
    allSub(idx) = '';
    hrate(idx) = [];
end

hitMean = mean(allHit)
hitSE   = std( allHit ) / sqrt( length( allHit ))

%% averaged over sessions
%  the previous data was on a single session level. because this is
%  difficult to put into a table, I now want to display these descriptors
%  on a subject level (including variance between sessions)
counter = 1;

skipit = zeros(size(supps,1),1);

pruned = {};
for ii = 1:size(supps,1)-1
    if skipit(ii) == 1
        continue
    end
    
    idx = strcmp(supps{ii,1}, supps(:,1));  % index of all sessions that belong to one subject
    pruned{counter,1}  = supps{ii,1};          % bidsID
    pruned{counter,2}  = sum(idx);             % number of sessions
    pruned{counter,3}  = mean([supps{idx,3}]); % NUMBER OF BUNDLES
    pruned{counter,4}  = mean([supps{idx,4}]); % NUMBER OF BUNDLES IN HIPP
    pruned{counter,5}  = mean([supps{idx,5}]);     % TRIAL NUMBER 5
    pruned{counter,6}  = std([supps{idx,5}]) / sqrt(length([supps{idx,5}]));
    
    pruned{counter,7}  = mean([supps{idx,6}]); % AND HITS 6
    pruned{counter,8}  = std([supps{idx,6}]) / sqrt(length([supps{idx,6}]));
    
    pruned{counter,9}  = mean([supps{idx,7}]); % NUMBER OF BUNDLES WITH SU 7
    pruned{counter,10} = std([supps{idx,7}]) / sqrt(length([supps{idx,7}]));
    
    pruned{counter,11} = mean([supps{idx,8}]); % NUMBER OF SU 8
    pruned{counter,12} = std([supps{idx,8}]) / sqrt(length([supps{idx,8}]));
    
    pruned{counter,13} = mean([supps{idx,9}]); % NUMBER OF PUTATIVE INDEX UNITS 9
    pruned{counter,14} = std([supps{idx,9}]) / sqrt(length([supps{idx,9}]));
    
    pruned{counter,15} = mean([supps{idx,10}]); % AND INDEX UNITS 10
    pruned{counter,16} = std([supps{idx,10}]) / sqrt(length([supps{idx,10}]));
    
    counter = counter +1;
    skipit(idx) = 1;
end

