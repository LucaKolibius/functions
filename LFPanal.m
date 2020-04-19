% bootstrapping to test for significance
% select MW with highest spike-field coupling
% throw out artefacts

%% lets start with each SU
IUpath = 'X:\Luca\indexSUvisu\varWindow_ENC_-1-113_RET_-1113_DP_p99_th1';
IUdir  = dir([IUpath, filesep, 'SU_*.mat']);
LFPdir = dir('X:\George\Analysis\Data Continuous\sub-*_micro_RAW_DS-1000_SPK-INT.mat');
artDir = dir('X:\Luca\lfpAnal\lfp_artefacts\LFPartefacts_*.mat');
% nonLFP_all = [];
% indLFP_all = [];

for it = 1:size(IUdir,1)
    clearvars -except IUpath IUdir LFPdir artDir it nonLFP_all indLFP_all
% it = 1;
load([IUdir(it).folder, filesep, IUdir(it).name]) % load IU info
load([artDir(1).folder, filesep, 'LFPartefacts_', bidsID, '_', sesh, '_enc_', wireName(1:end-1),'.mat'], 'artefactsEnc'); % load artefact timestamps (1x8 cell with each cell corresponding to a wire)
wireU = str2double(wireName(end));

% get old naming if you only have bids
if ~exist('subjID','var')
    subj = sub_ID_conversion(bidsID, 'yes');
    subjID = [subj, '_',sesh];
end

% get mSubject and mSession
mSubject = subjID(1:end-3);
mSession = subjID(end-1:end);
% if the session name is called 1b then this line prevents an error during cd
mSubject(regexp(mSubject,'_')) = [];
if isempty(regexp(mSession,'S', 'ONCE'))
    mSession = ['S', mSession];
end

% load encoding trigger
cd(['X:\Luca\data\', mSubject,filesep,  mSession,filesep]);
abc = dir('2*'); cd(abc.name);
p2d = [cd, filesep];
trls = 1; % can be ignored here; for when you want to have a specific animal cue
[~, hitsIdx, ~, ~, ~, ~, encTrigger, ~, ~] = loadLogs(p2d, trls);

% encoding trial +2s before and after
encTrigger = encTrigger(hitsIdx,1); % only hits, only cue-locked
encTrigger = round(encTrigger*1000); % seconds to samples

% load wirename
loadFile = dir(['X:\Luca\data\allSbj\allSU_', bidsID,'_',  mSession '.mat']);
load([loadFile.folder, filesep, loadFile.name])

idxWire = allSU{su,1};
idxBundle = idxWire(1:end-1);
disp(idxWire);

loadLFP = cellfun(@(x)regexp(x,[bidsID, '_', mSession]), {LFPdir.name}, 'UniformOutput', false);
loadLFP = find(~cellfun(@isempty, loadLFP));
load([LFPdir(loadLFP).folder, filesep, LFPdir(loadLFP).name])

% indexing the wires that I want to load
idx = cellfun(@(x)regexp(x,idxBundle), data_micro.label, 'UniformOutput', false);
idx = cellfun(@(x)find(x==1), idx, 'UniformOutput', false);
idx = ~cellfun(@isempty, idx);
idx = find(idx == 1);

data = [];
data.lfp_raw = data_micro.trial{1}(idx,:);
data.label = data_micro.label(idx);
data.idxTrl = sensTrlsFFPP;
data.tw = [encTrigger-1000 encTrigger];

% get artefacts out of raw (this is super ugly code)
for artMW = 1:size(data.lfp_raw,1) % for each MW
    artVec = []; % generate a vector with all the timestamps at which we have an artefact
    for t = 1:size(artefactsEnc{artMW},1)
        artVec = [artVec artefactsEnc{artMW}(t,1) : artefactsEnc{artMW}(t,2)];
    end
    tempVec = data.lfp_raw(artMW,:);
    tempVec(artVec) = NaN; % make these timepoints NaNs
    data.lfp_noArt(artMW,:) = tempVec;
end

for trl = 1:size(encTrigger,1)
    segmDat = data.lfp_noArt(:,data.tw(trl, 1):data.tw(trl, 2));
    delTrl = sum(isnan(segmDat),2) > 500; % more than half of the segment is an artefact
    segmDat(delTrl,:) = NaN; % chuck the whole segment
    data.segmDat{trl} = segmDat;
end

data.pnts  = size(data.segmDat{1},2);
data.srate = 1000;
data.chan  = size(data.label,1); 
    
%     cellfun(@(x)ismember(x, TRLtw), artefactsEnc, 'UniformOutput', false);
%     cellfun(@(x)ismember(TRLtw, x), artefactsEnc, 'UniformOutput', false);

% hz = linspace(0,1000,1001);
hz = linspace(0, data.srate, data.pnts);
nonLFP = []; indLFP = [];
% for mw = 1:8 % repeat for all 8 MW
for mw = wireU
    nonCounter = 0; % start at trial 1
    indCounter = 0;
    for trl = 1:size(encTrigger,1) % go over all trials
        if data.idxTrl(trl) == 0 % not indexed
            nonCounter = nonCounter + 1;
            %             nonLFP(mw, nonCounter,:) = abs(fft(data.segmDat{trl}(mw,:))/2001);
            nonLFP(mw, nonCounter,:) = abs(fft(data.segmDat{trl}(mw,:)) / data.pnts);
            
        elseif data.idxTrl(trl) == 1 % indexed trial
            indCounter = indCounter + 1;
%             indLFP(mw, indCounter,:) = abs(fft(data.segmDat{trl}(mw,:))/2001);
            indLFP(mw, indCounter,:) = abs(fft(data.segmDat{trl}(mw,:)) / data.pnts);
        end
    end
end

% % average over trials
% nonLFP = squeeze(nanmean(nonLFP,2));
% indLFP = squeeze(nanmean(indLFP,2));

nonLFP = squeeze(nonLFP(wireU,:,:));
indLFP = squeeze(indLFP(wireU,:,:));

% nonLFP_all(it,:) = nonLFP;
% indLFP_all(it,:) = indLFP;
nonLFP_all(it) = {nonLFP};
indLFP_all(it) = {indLFP};
end

%% find empirical difference
% find power spectrum difference per MW and then average
for mw = 1:size(IUdir,1)
    powDiff(mw,:) = nanmean(indLFP_all{mw},1) - nanmean(nonLFP_all{mw}, 1);
end

% create a variable that  has all the indLFP and nonLFP powerspectra in one
% double, separate for each MW
permLFP = {};
for i = 1:size(IUdir,1)
    permLFP(i) = {[indLFP_all{i}; nonLFP_all{i}]};
end

nperm = 10000;
for perm = 1:nperm % repeat 10.000 times
    for mw = 1:size(IUdir,1) % cycle through each IU
        permDx = randperm(size(permLFP{mw},1));
        powDiff_perm(mw,:,perm) = nanmean(permLFP{mw}(permDx(1:end-2),:),1) - nanmean(permLFP{mw}(permDx(end-1:end),:), 1);
    end
end

powDiff_perm = squeeze(nanmean(powDiff_perm,1));
diff_thUP = prctile(powDiff_perm,95,2);
diff_thLO = prctile(powDiff_perm, 5,2);
powDiff = nanmean(powDiff,1);

figure; hold on;
plot(powDiff, 'linew', 3)
xlim([0 100])
plot(diff_thUP, 'r--', 'linew', 1)
plot(diff_thLO, 'r--', 'linew', 1)
xlabel('Frequency [HZ]')
ylabel({'Power Difference Between'; 'Indexed and non-Indexed Trials'})
mAx = gca;
mAx.FontWeight = 'bold';
mAx.FontSize = 18;


for sub = 1:size(indLFP,1)
    subplot(4,2,sub)
    plot(hz, nonLFP(sub,:));
    hold on
    plot(hz, indLFP(sub,:));
    xlim([0 50])
end