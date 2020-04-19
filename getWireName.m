% getWireName
% gets the wire name on which the IU are
clear
datapathIU = 'X:\Luca\indexSUvisu\varWindow_ENC_0-213_RET_-1113_DP_p99_th1'; % IU
datapathSU = 'X:\Luca\data\allSbj'; % SU
numIU = dir([datapathIU, filesep, 'SU_*.mat']);
load('X:\Luca\indexSUvisu\varWindow_ENC_0-213_RET_-1113_DP_p99_th1\numSU.mat')

myVar = [];
allIU = [];
for ia = 1:size(numIU,1)
    % load in index neurons
    load([numIU(ia).folder, filesep, numIU(ia).name], 'bidsID', 'sesh', 'su')
    
    numSUdir = dir([datapathSU, filesep, 'allSU_sub-*.mat']);    

    lufor = [bidsID, '_', sesh]; % look for this sessio in numSU
    idx = ~cellfun(@isempty, cellfun(@(x) regexp(x,lufor), {numSUdir(:).name}, 'UniformOutput', false));
    
    load([numSUdir(idx).folder, filesep, numSUdir(idx).name])
    
myVar = [myVar, [bidsID; sesh; su; allSU(su,1)]];
allIU = [allIU; bidsID, sesh, allSU(su,:)];
end

