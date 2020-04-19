function fndLFP(ic)
%%
SUdatapath = 'X:\Luca\indexSUvisu\varWindow_ENC_-1-113_RET_-1113_DP_p99_th1\';
allIdxU = dir([SUdatapath, 'SU_*.mat']);

% [4 6:13] % takes overlap into account
% 1:size(allIdxU,1) % sometimes there is an overlap
% for ic = 1:size(allIdxU,1) % sometimes there is an overlap
load([allIdxU(ic).folder, filesep, allIdxU(ic).name], 'subjID', 'su', 'bidsID', 'sesh');

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
trls = 1; % can be ignored here
[~, hitsIdx, ~, ~, ~, retTrigger, encTrigger, ~, ~] = loadLogs(p2d, trls);

% encoding trial
encTrigger = encTrigger(hitsIdx,:); % only hits
encTrigger = round(encTrigger*1000); % seconds to samples

% % repeat for retrieval
% retTrigger = retTrigger(hitsIdx,:);
% retTrigger = round(retTrigger*1000);

% load wirename
loadFile = dir(['X:\Luca\data\allSbj\allSU_', bidsID,'_',  mSession '.mat']);
load([loadFile.folder, filesep, loadFile.name])

idxWire = allSU{su,1};
idxBundle = idxWire(1:end-1);
disp(idxWire);

% get relevant wires
% load in LFP data
% already downsampled to 1000 hz and spike interpolated
datapath = 'X:\George\Analysis\Data Continuous\';
dataFiles = dir([datapath, '*micro_RAW_DS-1000_SPK-INT.mat']);

%     if ~exist('bidsID', 'var'); bidsID = sub_ID_conversion(subjID, 'yes'); end

loadFile = cellfun(@(x)regexp(x,[bidsID, '_', mSession]), {dataFiles.name}, 'UniformOutput', false);
loadFile = find(~cellfun(@isempty, loadFile));
load([dataFiles(loadFile).folder, filesep, dataFiles(loadFile).name])

% indexing the wires that I want to load
idx = cellfun(@(x)regexp(x,idxBundle), data_micro.label, 'UniformOutput', false);
idx = cellfun(@(x)find(x==1), idx, 'UniformOutput', false);
idx = ~cellfun(@isempty, idx);
idx = find(idx == 1);

%% detect artefacts (manually)
% encoding
% filter linenoise

cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51; 49*2 51*2; 49*3 51*3];
cfg.bsfilttype = 'but'; % for erlangen
%     cfg.bsfilttype = 'fir';  % for bham
op_t = ft_preprocessing(cfg, data_micro); % filter linenoise


cfg = [];
cfg.viewmode = 'vertical';
cfg.trl = [encTrigger(:,1)-3*1000, encTrigger(:,3)+3*1000, zeros(size(encTrigger, 1))];
[encTrials] = ft_redefinetrial(cfg, op_t);
% [delTS_enc] = LFP_artDet(encTrials, 'X:\Luca\lfpAnal\lfp_autoAR', 8);
delTS_enc = encTrials;
% artefact rejection
artefactsEnc = [];
disp('enc');
for mw = 1:size(idx,1)
    disp(mw);
    
    % only select the current mircorwire
    cfg = [];
    cfg.channel = idx(mw);
    op_t = ft_selectdata(cfg, encTrials);
    
    cfg = [];
    %     cfg.artfctdef.visual.artifact = delTS_enc{idx(mw)};
        cfg.artfctdef.visual.artifact = artefactsEnc{mw};
    cfg.ylim = [-1200 1200];
    artRej = ft_databrowser(cfg, op_t);
    artefactsEnc{mw} = artRej.artfctdef.visual.artifact;
%     save(['X:\Luca\lfpAnal\lfp_artefacts\LFPartefacts_',bidsID,'_', sesh,'_enc_',idxBundle, '.mat'], 'artefactsEnc')
end
end
% X:\Luca\lfpAnal\lfp_artefacts
%     % retrieval
%     cfg = [];
%     cfg.viewmode = 'vertical';
%     cfg.trl = [retTrigger(:,1)-3*1000, retTrigger(:,3)+3*1000, retTrigger(:,1)];
%     [retTrials] = ft_redefinetrial(cfg, data_micro);
%         [delTS_ret] = LFP_artDet(retTrials, 'X:\Luca\lfpAnal\lfp_autoAR', 8);
%
%     % artefact rejection
%     artefactsRet = [];
%     disp('ret');
%     for mw = 1:size(idx,1)
%          disp(mw);
%         cfg = [];
%         cfg.channel = idx(mw);
%         op_t = ft_selectdata(cfg, retTrials);
%
%         cfg = [];
%         cfg.bsfilter = 'yes';
%         cfg.bsfreq = [49 51; 49*2 51*2; 49*3 51*3];
%         cfg.bsfilttype = 'but'; % for erlangen
%         %     cfg.bsfilttype = 'fir';  % for bham
%         op_t = ft_preprocessing(cfg, op_t); % filter linenoise
%
%         cfg = [];
%         cfg.artfctdef.visual.artifact = delTS_ret{idx(mw)};
%         cfg.ylim = [-600 600];
%         artRej = ft_databrowser(cfg, op_t);
%         artefactsRet{mw} = artRej.artfctdef.visual.artifact;
%     end
%
%     save(['X:\Luca\lfpAnal\lfp_artefacts\LFPartefacts_',bidsID,'_', sesh,'_ret_',idxBundle, '.mat'], 'artefactsRet')
% end
% save artefact cleansed files

% Gram-Schmidt orthogonalization of MW

% find MW with highest spike-field coupling (loop through bands)

% repeat process with un-orthogolized data if we do not find anything
% (indicative of theta playing a role)

% look if first spike is in a specific phase of a specific band
% look if first burst is in a specific phase of a specific band

% contrast this with non-indexed trials