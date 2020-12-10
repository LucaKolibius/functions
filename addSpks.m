clear
% allSbj = dir('X:\Luca\data\P*');
allSbj = dir('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\P*')
allSbj([1 2 4 5 11 ]) = [];
allSpk = struct;
counter = 0;
for sbj = 1 : length(allSbj)
    
    cd([allSbj(sbj).folder, filesep, allSbj(sbj).name])
    allSesh = dir('S*');
    
    for sesh = 1 : length(allSesh)
        
        if strcmp(allSbj(sbj).name, 'P05') && strcmp(allSesh(sesh).name, 'S5')
            continue % I couldn't split S5 from S4
        end
        
        if strcmp(allSbj(sbj).name, 'P13') && strcmp(allSesh(sesh).name, 'S4')
            continue % cannot yet get the trigger to match up
        end
        
        %% CD TO SESSION FOLDER
        cd([allSesh(sesh).folder, filesep, allSesh(sesh).name]);
        abc = dir('20*'); cd(abc.name);
        
        %% FIND SAMPLING RATE
        SR = findSR;
        
        %% LOAD LOGS
        [~, hitsIdx, ~, ~, ~, retTrigger, encTrigger, ~, ~] = loadLogs([cd, filesep], 1);
                
        %% FIND CLUSTER IN THE HIPPOCAMPUS
        cd advancedAnalysis\elecLoc
        abc = dir('elecLoc*'); load(abc.name);
        elecLoc(:,1) = regexprep(elecLoc(:,1), 'o_', '');
        hipBund = elecLoc(strcmp(elecLoc(:,2),'hipp'),1); % bundles in the hippocampus
        
        cd ..\manualRej
        abc = dir('clNames*'); load(abc.name);
        
        delCl = zeros(1,size(clNames,1));
        for cl = 1: size(clNames,1)
            
            whereis = [];
            for loc = 1: size(elecLoc,1)
                whereis(loc) = contains(clNames{cl}, elecLoc(loc,1)); % on which bundle in elecLoc is the current clName?
            end
            
            if ~any(whereis) % NO MATCH FOUND
                error('COULD NOT MATCH CLUSTER TO ANY BUNDLE POSITION!')
            end
            
            if ~strcmp(elecLoc(logical(whereis),2), 'hipp') % cluster is not on the hippocampus
                delCl(cl) = 1;
            end
            
        end
        clNames(logical(delCl)) = [];
               
        cd ..\postXCrej
        abc = dir('tableTimestamps*'); load(abc.name);
        p2d = cd;
        for cl = 1 : size(clNames,1)
            cd(p2d);
            isPos = contains(clNames(cl), 'Pos');
            clNum = str2double(clNames{cl}(end));
            
            if clNum == 0
                error('There must be 10 cluster on one wire')
            end
            
            wirename = clNames{cl}(1:end-6);
            wirename = regexprep(wirename,'times_Micro_', '');
            
            bundlename = wirename(1:end-1);
            
            spks = tableTimestamps_postXC.(clNames{cl});
            spks = spks{1};
            spks = spks/SR*1000; % SPIKES IN 1000HZ 
            
            bidsID = sub_ID_conversion(allSbj(sbj).name, 1);
            
            %% KICK OUT SOME ARTEFACTS
            if strcmp(bidsID, 'sub-0012') && strcmp(allSesh(sesh).name, 'S1') && strcmp(wirename, 'postHippL3') && length(spks) == 18
                continue % it's an artefact
            end
            
            if strcmp(bidsID, 'sub-1009') && strcmp(allSesh(sesh).name, 'S4') && strcmp(wirename, 'MB1') && length(spks) == 1184
                continue % it's an artefact
            end
            
            counter = counter + 1;
            allSpk(counter).bidsID      = bidsID;
            allSpk(counter).sesh        = allSesh(sesh).name;
            allSpk(counter).wirename    = wirename;
            allSpk(counter).bundlename  = bundlename;
            allSpk(counter).su          = clNum;
            allSpk(counter).isPos       = isPos;
            allSpk(counter).encTrigger  = encTrigger;
            allSpk(counter).retTrigger  = retTrigger;
            allSpk(counter).hitsIdx     = hitsIdx;
            allSpk(counter).spks        = spks;
            allSpk(counter).iu          = [];
            allSpk(counter).idxTrl      = [];
            allSpk(counter).initSR      = SR;
            
        end % END OF CLUSTER    
    end % END OF SESSION
end % END OF SBJ

allSpks = allSpk;
load('X:\Luca\data\allSbj\addtheseSpks.mat', 'addedSpks')

f = fieldnames(allSpks);
for i = 1:length(f)
    for entry = 1 : size(addedSpks,2)   
        allSpks(666+entry).(f{i}) = addedSpks(entry).(f{i});
    end
end
 
 save('X:\Luca\data\allSbj\allSpksHZ', 'allSpks')


%%
clear


% load('X:\Luca\data\allSbj\old\allSpks.mat', 'allSpks');
% oldSpks = allSpks; clear allSpks;
% 
% load('X:\Luca\data\allSbj\old\allSpksHZ_phs.mat', 'allSpks');
% phsSpks = allSpks; clear allSpks;
load('X:\Luca\data\allSbj\diffSpks.mat');
load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks');
newSpks = allSpks; clear allSpks;

%% Which spike was deleted in phsSpk?
% allPhase = {};
% for spk = 1: size(phsSpks,2)
%     allPhase(spk,1) = {[phsSpks(spk).bidsID, phsSpks(spk).sesh, phsSpks(spk).wirename, num2str(length(phsSpks(spk).spks))]};
% end

allNew = {};
for spk = 1: size(newSpks,2)
    allNew(spk,1) = {[newSpks(spk).bidsID, newSpks(spk).sesh, newSpks(spk).wirename, '_',num2str(length(newSpks(spk).spks))]};
end

allOld = {};
for spk = 1: size(oldSpks,2)
    allOld(spk,1) = {[oldSpks(spk).bidsID, oldSpks(spk).sesh, oldSpks(spk).wirename, '_',num2str(length(oldSpks(spk).spks))]};
end

%% % OLD vs PHS: 0012_S1_postHippL3 17 needs to go(?)



%% In old variable, but missing from new variable
isthere = [];
for spk = 1 : length(allOld)
    isthere(spk) = any(strcmp(allOld(spk), allNew));
end
tt= find(isthere == 0)
allOld(tt)

allOld = {};
for spk = 1: size(oldSpks,2)
    allOld(spk,1) = {[oldSpks(spk).bidsID, oldSpks(spk).sesh, oldSpks(spk).wirename, '_', num2str(oldSpks(spk).su), '_spks', num2str(length(oldSpks(spk).spks))]};
end

missing = allOld(tt)

%% adding the missing spikes
addpath(genpath('X:\Luca\toolboxes\wave_clus-master')); % waveclus 3.0
add = 1 : size(missing,1)
    temp = missing{add};
    bidsID = temp(1:8);
    sesh   = temp(9:10); sesh = regexprep(sesh, '_S', 'S1b');
    wire = temp(11:regexp(temp,'_')-1);
    su = temp(regexp(temp,'_')+1);
    spks = temp(regexp(temp, '_spks')+5:end);
    subjID = sub_ID_conversion(bidsID,1);
    
    cd(['X:\Luca\data\', subjID, filesep, sesh])
    abc = dir('20*'); cd(abc.name);
    
%     switch isPos
%         case 0
%             cd('..\negDetect')
%         case 1
            cd('posDetect')
%     end
    
    pnClus = dir(['times_*', wire,'.mat']); % list of all electrode names
    wave_clus(pnClus.name); % loads results from automatic clustering

        %% In new variable, but missing from old one (like a few new negative cluster)
        %% lets see who needs to go
        %  (the ones that were falsely added!)
        
        isthere = [];
        for spk = 1 : length(allNew)
            isthere(spk) = contains(allNew(spk), allOld);
        end
        tt = find(isthere == 0);
        
        allNew = {};
        for spk = 1: size(newSpks,2)
            allNew(spk,1) = {[newSpks(spk).bidsID, newSpks(spk).sesh, newSpks(spk).wirename, '_', num2str(newSpks(spk).su), '_isPos', num2str(newSpks(spk).isPos), '_spks', num2str(length(newSpks(spk).spks))]};
        end
        tooMuch = allNew(tt);
        
        if length(tooMuch) ~= 115
            error('loading problem');
        end
    tooMuch(55:114) = [];
    tooMuch([8 10 41:44 50 51 53 54 55]) = []; % 8 is a new negative cluster, 10 is a fine positive cluster, 41 are the 1bs, 44 could be an artefact (hatty), 
    % 50 and 51 are normal cluster, 55 is a 1b, 53 and 54 are negative
    % cluster
    addpath(genpath('X:\Luca\toolboxes\wave_clus-master')); % waveclus 3.0
%     rmv = 1 : size(tooMuch,1)
    
    rmv = rmv + 1
    temp = tooMuch{rmv};
    bidsID = temp(1:8);
    sesh   = temp(9:10); sesh = regexprep(sesh, '_S', 'S1b');
    wire = temp(11:regexp(temp,'_')-1); wire = regexprep(wire, 'b', '');
%     su = temp(regexp(temp,'_')+1);
    spks = temp(regexp(temp, '_spks')+5:end);
    subjID = sub_ID_conversion(bidsID,1);
    isPos = str2double(temp(regexp(temp, 'isPos')+5));
    cd(['X:\Luca\data\', subjID, filesep, sesh])
    abc = dir('20*'); cd(abc.name);
    
    switch isPos
        case 0
            cd('negDetect')
        case 1
            cd('posDetect')
    end
    
    pnClus = dir(['times_*', wire,'.mat']); % list of all electrode names
    wave_clus(pnClus.name); % loads results from automatic clustering
    
    
        pnClus = dir(['times_*', 'MRTA2','.mat']); % list of all electrode names
    wave_clus(pnClus.name); % loads results from automatic clustering
