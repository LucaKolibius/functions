% handles the load-in of data for visuSU
function [hitsIdx, subj, spksEnc, spksRet, sensTrlsFFPP, animalCues, WS, averageWS, ffTrls, ppTrls, ewp, associates] = loadInDat(allSpks, su)
sensTrlsFFPP = [];
ffTrls = [];
ppTrls = [];

global prePath;
cd([prePath, '\Luca\data']);

bidsID          = allSpks(su).bidsID;
subj = sub_ID_conversion(bidsID, 'yes');

sesh            = allSpks(su).sesh;
sesh = regexprep(sesh, 'S1b', 'S1');
sr              = allSpks(su).initSR;

encTrig         = allSpks(su).encTrigger(:,[1:3]); % in seconds
retTrig         = allSpks(su).retTrigger(:,[1:3]); % in seconds
retRT           = allSpks(su).retRT;
idxTrl          = allSpks(su).idxTrl;
hitsIdx         = allSpks(su).hitsIdx;
wirename        = allSpks(su).wirename;

cd(subj)
cd(sesh)
abc = dir; cd(abc(3).name);

% ANIMAL CUES
[~, ~, ~, ~, ~, ~, ~, ~, ~, animalCues, associates] = loadLogs([cd, filesep], idxTrl);
% load single units from processed variable
% load single units from raw data (old)
clusterSpikes = allSpks(su).spks/1000;

% adjusted for ENC(-1-113) RET(-1113)
spksEnc = insertSpiketimes2(encTrig, clusterSpikes, [1 3], [-3 0])'; % 3 seconds prior to cue trigger until end of trial
spksRet = insertSpiketimes2(retTrig, clusterSpikes, [1 3], [-3 0])'; % 3 seconds prior to cue trigger until end of trial

% This is to visualize the reinstatement value that goes from enc (cue+2 to
% end) and ret (cue to response)
[ewp.enc, ewp.twEnc] = insertSpiketimes2([encTrig(:,1) encTrig(:,3)], clusterSpikes, [1 2], [0 0]);
[ewp.ret, ewp.twRet] = insertSpiketimes2([retTrig(:,1)   retRT(:,1)],   clusterSpikes, [1 2], [0 0]);

%% WAVESHAPE
switch allSpks(su).isPos
    case 0
        posNeg = 'neg';
    case 1
        posNeg = 'pos';
end

if exist('detect_newSR', 'file') == 7
    cd('detect_newSR');
else
    cd([posNeg, 'Detect']); % go to posDetect or negDetect
end

% LOAD WAVESHAPE FROM DETECT FOLDER OR detect_newSR
abc = dir(['times_*', wirename, '.mat']);
if size(abc,1)>1; error('Naming conflict before loading file'); end
load(abc.name, 'spikes', 'cluster_class', 'par')

clIdx        = cluster_class(:,1) == allSpks(su).su;
numSpksVar   = sum(clIdx);
numSpks      = size(clusterSpikes,1);

%% IF SPIKETIMES ARE NOT CORRECT IT IS PROBABLY IN detect_newSR
if ~isequal(par.sr, sr) | ~isequal(numSpksVar, numSpks)
    cd('..\detect_newSR2');
    abc = dir(['times_*', wirename,'_' posNeg, '.mat']);
    if size(abc,1)>1; error('Naming conflict before loading file'); end
    load(abc.name, 'spikes', 'cluster_class', 'par', 'inspk')
    
    % re-calculate number of spikes
    clIdx        = cluster_class(:,1) == allSpks(su).su;
    numSpksVar   = sum(clIdx);
end

if ~isequal(par.sr, sr) | ~isequal(numSpksVar, numSpks)
    error('Still not load the correct waveclus file')
end

WS = spikes(clIdx,:);
averageWS = mean(WS,1);


%     %% load single units from raw data (old)
%     cd advancedAnalysis\postXCrej
%     loadVar = dir('tableTimestamps_postXC_*'); load(loadVar.name); % load tableTimestamps_postXC
%     cd ../manualRej/
%     loadVar = dir('clNames*'); load(loadVar.name); % load clNames
%     cd ../..
%     clCount = 0;
%     allSpikes = [];
%     for i=1:size(tableTimestamps_postXC,2)
%         mVarname = tableTimestamps_postXC.Properties.VariableNames{1,i};
%         if ismember(mVarname, clNames)
%             clCount = clCount +1;
%             clusterSpikes = table2array(tableTimestamps_postXC{1,i})'; % load in all spiketimes from cluster mVarname
%             clusterSpikes = clusterSpikes/32000; % samples to sec
%
%             % segment trials of this SU into trials
%             % old: here I inserted spikes in a time window surrounding one
%             % trigger
%             %             spksEnc(:, clCount) = insertSpiketimes(encTrigger, clusterSpikes, 1, [-2 5])'; % cue locked; 6s  %I could also generate a big file with all subj data in it (using a struct) or put enc and ret in one variable, but what for?
%             %             spksRet(:, clCount) = insertSpiketimes(retTrigger, clusterSpikes, 1, [-2 3])'; % cue locked; 4s
%
%             % adjusted for ENC(0-213) RET(-1113)
%             spksEnc(:, clCount) = insertSpiketimes2(encTrigger, clusterSpikes, [1 3], [-2 0])'; % 2 seconds prior to cue trigger until response trigger
%             spksRet(:, clCount) = insertSpiketimes2(retTrigger, clusterSpikes, [1 3], [-3 3])'; % 3 seconds prior to cue trigger until 3 seconds post response trigger
%
%             % get waveshape
%             demGate = cd; % where the ncs files are
%
%             if ~isempty(regexp(mVarname, 'Pos', 'ONCE'))
%                 posNeg = 'posDetect';
%             elseif ~isempty(regexp(mVarname, 'Neg', 'ONCE'))
%                 posNeg = 'negDetect';
%             end
%
%             cd(posNeg)
%
%             temp1 = mVarname;
%             temp1(end-5:end) = [];
%
%             % fixes a problem with the file names yet again.
%             if size(temp1,2) == 6
%                 temp1(1:2) = [];
%             end
%             temp = dir(['times_CSC_', temp1, '.mat']);
%
%             % problem due to different file names. this fixes it. Its not
%             % pretty, but does the job.
%             if isempty(temp)
%                 temp = dir(['times_Micro_', temp1, '.mat']);
%             end
%
%             load(temp.name, 'spikes', 'cluster_class');
%
%             clNum = mVarname(end);
%             extrCl = cluster_class(:,1)==str2num(clNum);
%
%             counterX = 0;
%             for ib = 1 : size(extrCl,1)
%                 if extrCl(ib,1) == 1
%                     counterX = counterX+1;
%                     allSpikes(clCount).waveshape(counterX, :) = spikes(ib,:);
%                 end
%             end
%             averageWS(clCount,:) = mean(allSpikes(clCount).waveshape,1);
%             cd(demGate)
%         end
%     end
end % END OF FUNCTION