%% This function gets the absolute number of ff/fp/pp hits & ff/fp/pp misses
% the difference here is that numHit2 and numMiss2 does not add the previous values for each new value
function [numHit2, numMiss2] = get_FfFpPp_hits_miss(subjID)
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
cd advancedAnalysis\finalDat\spikenumbers\
abc = dir('encSpikeNumber_cueLocked*');
load(abc.name);

% hit
numFF_hit = 0;
numPP_hit = 0;
numFP_hit = 0;
for i=1:size(encSpikeNumber_cueLocked,2)
    if strcmp(encSpikeNumber_cueLocked{1,i}, 'ff') && strcmp(encSpikeNumber_cueLocked{2,i}, 'hit')
        numFF_hit=numFF_hit+1;
    
    elseif strcmp(encSpikeNumber_cueLocked{1,i}, 'pp') && strcmp(encSpikeNumber_cueLocked{2,i}, 'hit')
        numPP_hit=numPP_hit+1;
    
    elseif strcmp(encSpikeNumber_cueLocked{1,i}, 'fp') && strcmp(encSpikeNumber_cueLocked{2,i}, 'hit')
        numFP_hit=numFP_hit+1;
    end
end

% miss
numFF_miss = 0;
numPP_miss = 0;
numFP_miss = 0;
for i=1:size(encSpikeNumber_cueLocked,2)
    if strcmp(encSpikeNumber_cueLocked{1,i}, 'ff') && strcmp(encSpikeNumber_cueLocked{2,i}, 'miss')
        numFF_miss=numFF_miss+1;
    
    elseif strcmp(encSpikeNumber_cueLocked{1,i}, 'pp') && strcmp(encSpikeNumber_cueLocked{2,i}, 'miss')
        numPP_miss=numPP_miss+1;
    
    elseif strcmp(encSpikeNumber_cueLocked{1,i}, 'fp') && strcmp(encSpikeNumber_cueLocked{2,i}, 'miss')
        numFP_miss=numFP_miss+1;
    end
end

numHit2  = [numFF_hit numPP_hit numFP_hit]; % the difference here is that numHit2 and numMiss2 does not add the previous values for each new value
numMiss2 = [numFF_miss numPP_miss numFP_miss];

end % end of function
