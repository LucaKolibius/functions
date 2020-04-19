function [encSpiketimes_cueLocked, encSpiketimes_respLocked, encSpiketimes_preTrial, encSpiketimes_stimLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID)
try
    cd X:/Luca/data
catch
    cd /media/ldk898/rds-share/Luca/data
end

mSubject = subjID(1:end-3);
mSession = subjID(end-1:end);

% if the session name is called 1b then this line prevents an error during cd
mSubject(regexp(mSubject,'_')) = []; 
if isempty(regexp(mSession,'S', 'ONCE'))
    mSession = ['S', mSession];
end

cd(mSubject)
cd(mSession)
abc = dir;
cd(abc(3).name)
p2d = cd;
p2d(end+1)='\';
[tableTemplate, ~, ~, ~, ~, retTrigger, encTrigger, ~, ~] = loadLogs(p2d, 1);

cd advancedAnalysis
cd postXCrej
abc1 = dir('tableTimestamps_postXC_*');
load(abc1.name)

% spiketimes tables
encSpiketimes_cueLocked = tableTemplate;
encSpiketimes_respLocked = tableTemplate;
retSpiketimes_cueLocked = tableTemplate;
retSpiketimes_respLocked = tableTemplate;

% spikenumber tables
encSpikeNumber_cueLocked  = tableTemplate;
encSpikeNumber_respLocked = tableTemplate;
retSpikeNumber_cueLocked  = tableTemplate;
retSpikeNumber_respLocked = tableTemplate;

%% add two more time windows for encoding
% preTrial
encSpikeNumber_preTrial    = tableTemplate;
encSpiketimes_preTrial     = tableTemplate;

% stimLocked
encSpikeNumber_stimLocked = tableTemplate;
encSpiketimes_stimLocked  = tableTemplate;

%% timewindows
timeWindow = [];
timeWindow.encCueLocked  = [-1 5];
timeWindow.encRespLocked = [-2 1];
timeWindow.retCueLocked  = [-1 2];
timeWindow.retRespLocked = [-2 1];
timeWindow.encPreTrial   = [-1 0]; % from cue
timeWindow.encStimLocked = [-1 3]; % from stim

clNames = {};
for iy = 1:size(tableTimestamps_postXC,2)
    varname = tableTimestamps_postXC.Properties.VariableNames{iy};
    
    if strcmp(varname(end), '0')
        continue
    end
    
    clNames(1, size(clNames,2)+1) = cellstr(varname);
    
    a = table2array(tableTimestamps_postXC{1,iy});
    a = a/32000; % samples to secs
    
    % add spiketimes of the cluster to the table
    % Encoding + CueLocked
    locking = 1; % for CueLocked
    encSpiketimes_cueLocked{size(encSpiketimes_cueLocked,1)+1,:}   = insertSpiketimes(encTrigger,a,locking,timeWindow.encCueLocked);
    
    % Encoding + RespLocked
    locking = 3; % for RespLocked
    encSpiketimes_respLocked{size(encSpiketimes_respLocked,1)+1,:} = insertSpiketimes(encTrigger,a,locking,timeWindow.encRespLocked);
    
    % Encoding + PreTrial
    locking = 1;
    encSpiketimes_preTrial{size(encSpiketimes_preTrial,1)+1,:}     = insertSpiketimes(encTrigger,a,locking,timeWindow.encPreTrial);
    
    % Encoding + StimLocked
    locking = 2;
    encSpiketimes_stimLocked{size(encSpiketimes_stimLocked,1)+1,:} = insertSpiketimes(encTrigger,a,locking,timeWindow.encStimLocked);
    
    % Retrieval + CueLocked
    locking = 1;
    retSpiketimes_cueLocked{size(retSpiketimes_cueLocked,1)+1,:}   = insertSpiketimes(retTrigger,a,locking,timeWindow.retCueLocked);
    
    % Retrieval + RespLocked
    locking = 3;
    retSpiketimes_respLocked{size(retSpiketimes_respLocked,1)+1,:} = insertSpiketimes(retTrigger,a,locking,timeWindow.retRespLocked);
    
end
clNames=clNames';
% beep
% open encSpiketimes_cueLocked
% open encSpiketimes_respLocked
% open encSpiketimes_preTrial
% open encSpiketimes_stimLocked
% 
% open retSpiketimes_cueLocked
% open retSpiketimes_respLocked
% open clNames
mCD = cd;
end