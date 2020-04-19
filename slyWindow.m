%% Sliding Window Approach
clear ENC
ENC.zehnMS(1,1).trial = zeros(5,5);
ENC.zehnMS(1,2).trial = zeros(5,5);
ENC.zehnMS(1,3).trial = zeros(5,5);

clear ENC
ENC.zehnMS(1,1).trial = zeros(5,5);
ENC.zehnMS(2,1).trial = zeros(5,5);
ENC.zehnMS(3,1).trial = zeros(5,5);
%%

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
[tableTemplate, ~, ~, ~, ~, retTrigger, encTrigger] = loadLogs(p2d);

cd advancedAnalysis
cd postXCrej
abc1 = dir('tableTimestamps_postXC_*');
load(abc1.name)

% spikenumber tables
encSpikeNumber_cueLocked  = tableTemplate;
retSpikeNumber_respLocked = tableTemplate;

%% 
WIP