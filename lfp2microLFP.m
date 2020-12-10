% % various settings
% LFPfold = 'X:\George\Analysis\Data Continuous\no-SPK-interpolation';                     % LFP;
% LFPdir  = dir('X:\George\Analysis\Data Continuous\sub-*.mat.mat');
% artDir = 'Z:\hanslmas-ieeg-compute\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS
% 
% % start with each SU
% for it = 1:size(LFPdir,1)
%     clearvars -except it LFPfold LFPdir artDir
%     disp(it)
%        
%     % get current bidsID & sesh
%     bidsID = LFPdir(it).name; sesh = bidsID;
%     bidsID = bidsID(1:8);
%     sesh   = sesh(10:11);
%         
%     load([artDir, bidsID, '_', sesh, '_micro_th-8.mat'], 'delIx', 'del_sampleinfo'); % load artefact timestamps (1x8 cell with each cell corresponding to a wire)
%     load([LFPfold, bidsID, '_', sesh, '_micro_RAW_DS-1000_SPK-INT.mat'], 'data_micro') 
%     
%     labs    = data_micro.label; % all labels
%     noSpks  = cellfun(@(x) regexp(x, 'SPKS'), labs, 'un', 0); % all labels apart from the spikes
%     noTrig  = cellfun(@(x) regexp(x, 'TRIG'), labs, 'un', 0); % indexes the trigger
%     noHits  = cellfun(@(x) regexp(x, 'enc+'), labs, 'un', 0); % indexes the hits / miss info
%     selChan = cellfun(@isempty, noSpks);
%     noTrig  = ~cellfun(@isempty, noTrig);
%     noHits  = ~cellfun(@isempty, noHits);
%     
%     % index all the names of LFPs (ignores spikes / trigger / hits and misses)
%     selChan      = logical(selChan - noTrig - noHits);   
%     cfg          = [];
%     cfg.channel  = data_micro.label(selChan);       % all 8 channels
%     data         = ft_selectdata(cfg, data_micro);  % select
%     
%     
%     for chan = 1 : size(del_sampleinfo,2)
%         artChan = del_sampleinfo{chan};
%         for art = 1 : size(artChan,1)
%             data.trial{1}(chan, artChan(art,1) : artChan(art,2) ) = nan;
%         end
%     end
%     
%     save(['X:\Luca\data\microLFP\', bidsID, '_', sesh, '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data')
%     
% end

%% FROM LAB PC
% various settings
% CHANGED FOR NO SPK INTERPOLATION!!

clear
addpath('X:\Common\toolboxes\fieldtrip-20200310');
ft_defaults;
% LFPfold = 'X:\George\Analysis\Data Continuous\no-SPK-interpolation\';   % LFP;
LFPfold = '\\analyse4.psy.gla.ac.uk\Project0309\Luca\Georges LFP\';
LFPdir  = dir([LFPfold, 'sub-*_micro_RAW_DS-1000.mat']);
% artDir = 'X:\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS


%% ONLY SELECT NEW DATA
newLFPidx = cell2mat(cellfun(@(x) contains(x, '01-Dec-2020'), {LFPdir.date}, 'un', 0));
LFPdir = LFPdir(newLFPidx);
% start with each SU
for it = 1:size(LFPdir,1)
    clearvars -except it LFPfold LFPdir artDir
    disp(it)
       
    % get current bidsID & sesh
    bidsID = LFPdir(it).name; sesh = bidsID;
    bidsID = bidsID(1:8);
    sesh   = sesh(10:11);
        
%     load([artDir, bidsID, '_', sesh, '_micro_th-8.mat'], 'delIx', 'del_sampleinfo'); % load artefact timestamps (1x8 cell with each cell corresponding to a wire)

%% MICRO
load([LFPfold, bidsID, '_', sesh, '_micro_RAW_DS-1000.mat'], 'data_micro') 

%% MACRO
%     load([LFPdir(it).folder, filesep, LFPdir(it).name], 'data_macro')
    
    %% FOR MICROS
    labs    = data_micro.label; % all labels
    noSpks  = cellfun(@(x) regexp(x, 'SPKS'), labs, 'un', 0); % all labels apart from the spikes
    noTrig  = cellfun(@(x) regexp(x, 'TRIG'), labs, 'un', 0); % indexes the trigger
    noHits  = cellfun(@(x) regexp(x, 'enc'), labs, 'un', 0); % indexes the hits / miss info
    selChan = cellfun(@isempty, noSpks);
    noTrig  = ~cellfun(@isempty, noTrig);
    noHits  = ~cellfun(@isempty, noHits);
    
%     % index all the names of LFPs (ignores spikes / trigger / hits and misses)
    selChan      = logical(selChan - noTrig - noHits);   
    
    %% FOR MACROS
%     labs = cellfun(@(x) x(end), data_macro.label, 'un', 0);
%     labs = str2double(labs);    
%     selChan = labs == 1;
    
    cfg          = [];
    cfg.channel  = data_micro.label(selChan);       % all 8 channels
    data         = ft_selectdata(cfg, data_micro);  % select
    
    
%     for chan = 1 : size(del_sampleinfo,2)
%         artChan = del_sampleinfo{chan};
%         for art = 1 : size(artChan,1)
%             data.trial{1}(chan, artChan(art,1) : artChan(art,2) ) = nan;
%         end
%     end
    
save(['\\analyse4.psy.gla.ac.uk\Project0309\Luca\data\microLFP\', bidsID, '_', sesh, '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data')     

end