%% YOU HAVE TO INCLUDE SBJ = 15 (P08ERL) 
clear all

allSbj = dir('Z:\hanslmas-ieeg-compute\Luca\data\P*');
load('X:\Luca\data\allSbj\allSpks.mat', 'allSpks')
load('X:\Luca\data\allSbj\rpplBund.mat', 'rpplBund')
counter = 1;

try
for sbj = 1 : length(allSbj)
    
    % these are empty. skip them
    if any(sbj == [1 2 4 5 11 15 17])
        continue
    end
    bidsID = sub_ID_conversion(allSbj(sbj).name, 1);
    
    cd([allSbj(sbj).folder, filesep, allSbj(sbj).name])
    allSesh = dir('S*');
    
    supps{counter, 1} = allSbj(sbj).name;
    
    for sesh = 1:length(allSesh)
        
        % skip these (no SU)
        if or( (sbj == 8 && sesh == 1), (sbj == 8 && sesh == 3) )
            continue
        end
        
        % skip this (cant untangle sesh 5 from sesh 4)
        if (sbj == 9 && sesh == 5)
            continue
        end
        
        supps{counter, 2} = allSesh(sesh).name;
         cd([allSesh(sesh).folder, filesep, allSesh(sesh).name])
            abc = dir('20*'); cd(abc(1).name);
            cd('advancedAnalysis\')
            cd('elecLoc')

        
        abc = dir('elecLoc*');
        load(abc.name)
        
        supps{counter, 3} = size(elecLoc,1); % number of bundles
        supps{counter, 4} = sum(strcmp(elecLoc(:,2), 'hipp')); % number of bundles in hipp
        
        
        sameSesh = and(contains({allSpks.bidsID}, bidsID), strcmp({allSpks.sesh}, allSesh(sesh).name) ); % same subject + session + bundle
        
        
        idx = find(sameSesh==1); idx = idx(1);
        supps{counter, 5} = size(allSpks(idx).encTrigger, 1);
        supps{counter, 6} = size(allSpks(idx).hitsIdx, 1);

        %% number of bundles with SU
        allBund = {allSpks(sameSesh).wirename};
        allBund = cellfun(@(x) x(1:end-1), allBund, 'un', 0);
        allBund = unique(allBund);
        
            supps{counter, 7} = length(allBund);

        %% number of SU
        supps{counter,8} = sum(sameSesh);
        
        %% number of putative index units
        iu = [allSpks(sameSesh).iu];
            supps{counter, 9} = sum(iu == 1);

        
        %% number of  index units
        iu = [allSpks(sameSesh).iu];
            supps{counter, 10} = sum(iu == 2);
        
        counter = counter + 1;
    end
end
catch e
     fprintf(2,'There was an error! The message was:\n%s \n',e.message);
        
    disp(sbj)
    disp(bidsID)
    disp(sesh)
end