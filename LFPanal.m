IUdir  = dir('X:\Luca\indexSUvisu\varWindow_ENC_-1-113_RET_-1113_DP_p99_th1\SU_*.mat';);
LFPdir = dir('X:\George\Analysis\Data Continuous\sub-*_micro_RAW_DS-1000_SPK-INT.mat');
artDir = dir('X:\Luca\engram allocation lfp\lfpArtefacts\lfp_artefacts\LFPartefacts_*.mat');
addpath('X:\Luca\toolboxes\gedCFC_tutorial');

% start with each SU
for it = 1%:size(IUdir,1)
    clearvars -except IUpath IUdir LFPdir artDir it nonLFP_all indLFP_all
    load([IUdir(it).folder, filesep, IUdir(it).name]) % load IU info
    load([artDir(1).folder, filesep, 'LFPartefacts_', bidsID, '_', sesh, '_enc_', wireName(1:end-1),'.mat'], 'artefactsEnc'); % load artefact timestamps (1x8 cell with each cell corresponding to a wire)
    wireU = str2double(wireName(end));
    
    % get old naming if you only have bids
    if ~exist('subjID','var')
        subj   = sub_ID_conversion(bidsID, 'yes');
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
    abc  = dir('2*'); cd(abc.name);
    p2d  = [cd, filesep];
    trls = 1; % can be ignored here; for when you want to have a specific animal cue
    [~, hitsIdx, ~, ~, ~, ~, encTrigger, ~, ~] = loadLogs(p2d, trls);
    
    % encoding trial
    encTrigger = encTrigger(hitsIdx,1); % only hits, only cue-locked
    encTrigger = round(encTrigger*1000); % seconds to samples
    
    % load wirename
    loadFile = dir(['X:\Luca\data\allSbj\allSU_', bidsID,'_',  mSession '.mat']);
    load([loadFile.folder, filesep, loadFile.name])
    
    idxWire   = allSU{su,1};
    idxBundle = idxWire(1:end-1);
    disp(idxWire);
    
    loadLFP = cellfun(@(x)regexp(x,[bidsID, '_', mSession]), {LFPdir.name}, 'UniformOutput', false); % find the LFP date for the correct patient and session ...
    loadLFP = find(~cellfun(@isempty, loadLFP));
    load([LFPdir(loadLFP).folder, filesep, LFPdir(loadLFP).name]) % ... and load it in
    
    idx = cellfun(@(x)regexp(x,idxBundle), data_micro.label, 'UniformOutput', false); % find the position of the MW of my bundle of interest
    idx = cellfun(@(x)find(x==1), idx, 'UniformOutput', false);
    idx = ~cellfun(@isempty, idx);
    idx = find(idx == 1);
    
    data = [];
    data.lfp_raw = data_micro.trial{1}(idx,:);
    data.label   = data_micro.label(idx);
    data.idxTrl  = sensTrlsFFPP;
    data.tw      = [encTrigger-1000 encTrigger];
    
    % reject artefacts
    for artMW  = 1 : size(data.lfp_raw,1) % for each MW
        artVec = []; % generate a vector with all the timestamps at which we have an artefact
        for t = 1:size(artefactsEnc{artMW},1) % create a vector with timepoints of all artefact sample points
            artVec = [artVec artefactsEnc{artMW}(t,1) : artefactsEnc{artMW}(t,2)];
        end
        tempLFP                 = data.lfp_raw(artMW,:);
        tempLFP(artVec)         = NaN; % make these timepoints NaNs
        data.lfp_noArt(artMW,:) = tempLFP; % save this LFP as artefact rejected LFK (lfp_noArt)
    end
    
    % next we look at trials and if they contain more than a predefined number
    % of nans. If they do, we chuck the whole trial
    for trl = 1:size(encTrigger,1)
        segmDat           = data.lfp_noArt(:,data.tw(trl, 1):data.tw(trl, 2));
        delTrl            = sum(isnan(segmDat),2) > 500; % more than half of the segment is an artefact
        segmDat(delTrl,:) = NaN; % chuck the whole segment
        data.segmDat{trl} = segmDat;
    end
    
    data.pnts  = size(data.segmDat{1},2);
    data.srate = 1000;
    data.chan  = size(data.label,1);
    
    % Get the spikes
    idxSpks = allSU{su,3}/32;
    
    nwin    = 50; % time in samples before and after spikes that is considered for analysis
    filtf   = 0; % highpass filter at 10hz; set to 0 for no filtering
    [evec, signSFC] = gedSFC_method1(data, idxSpks, nwin, filtf, bidsID, sesh, su);
    
    if signSFC == 0
        filtf   = 10; % highpass filter at 10hz; set to 0 for no filtering
        gedSFC_method1(data, idxSpks, nwin, filtf, bidsID, sesh, su) % redo SFC
    end
    
    % PPC_linenoise(data, idxSpks) % if I find some spikelocking I wont know if it is significant, add sine to my spikes and rerun method1?
    
    % beta peak because padding is too small
    % clearvars -except data idxSpks bidsID sesh su it % memory problems
    % % Do SFC for the continuous data segmented into 1s snippets (to combat NaNs) until the end
    % % of the recording
    % gedSFC_segm    (data, idxSpks, bidsID, sesh, su);
    %
    
    % % downsample first to 100hz
    % datadown.trial = {data.lfp_raw};
    % datadown.time = data_micro.time;
    % datadown.label = {'Chan1'; 'Chan2'; 'Chan3'; 'Chan4'; 'Chan5'; 'Chan6'; 'Chan7'; 'Chan8'};
    % cfg = [];
    % cfg.resamplefs = 300; % resample at 100hz
    % cfg.detrend = 'no';
    % cfg.trials = 'all';
    %
    % datadown = ft_resampledata(cfg, datadown);
    %
    % idxSpks = allSU{su,3}/32000*300;
    % inptDat.lfp_noArt = datadown.trial{1};
    % gedSFC_segm    (inptDat, idxSpks, bidsID, sesh, su);
    
    % % Do SFC for the trial segmeneted data
    % gedSFC_trlsegm (downdata, idxSpks, bidsID, sesh, su);
    
if signSFC == 1
data.segmComp = cellfun(@(x) x'*evec, data.segmDat, 'un', 0); 
numtrl = size(data.idxTrl,2);
nonLFP = []; indLFP = [];
   
nonCounter = 0; % start at trial 1
indCounter = 0;

for trl = 1:numtrl
    if data.idxTrl(trl) == 0 % not indexed
        nonCounter = nonCounter + 1;
        nonLFP(nonCounter,:) = abs(fft(data.segmComp{trl}) / data.pnts);
        
    elseif data.idxTrl(trl) == 1 % indexed trial
        indCounter = indCounter + 1;
        indLFP(indCounter,:) = abs(fft(data.segmComp{trl}) / data.pnts);
    end
end

% nonLFP_all(it) = {nonLFP};
% indLFP_all(it) = {indLFP};

%% find empirical difference
% find power spectrum difference per MW and then average
powDiff(it,:) = nanmean(nonLFP,1) - nanmean(indLFP, 1);

% create a variable that  has all the indLFP and nonLFP powerspectra in one
% double, separate for each MW

permLFP{it} = [nonLFP; indLFP];
numIdx(it) = sum(data.idxTrl);

savepath = ['X:\Luca\engram allocation lfp\SFC\method1\', num2str(filtf), filesep, num2str(nwin), filesep];
% savepath = ['/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/fileshare/engram allocation lfp/', num2str(filtf), '/', num2str(nwin), '/'];
save([savepath,  'LFP_', bidsID, '_', sesh, '_su', num2str(su)], 'numIdx', 'permLFP', 'powDiff')
end % if condition end (if I have significant SFC)

end

% nperm = 10000;
% for perm = 1:nperm % repeat 10.000 times
%     for iu = 1:size(IUdir,1) % cycle through each IU
%         
%         if isempty(permLFP{iu}) % no significant SFC
%             continue
%         end
%         
%         numIdxTrl               = numIdx{iu}; % number of trials this unit indexes
%         permDx                  = randperm(size(permLFP{iu},1));
%         powDiff_perm(iu,:,perm) = nanmean(permLFP{iu}(permDx(1:end-numIdxTrl),:),1) - nanmean(permLFP{iu}(permDx(end-numIdxTrl+1:end),:), 1);
%     end
% end
% 
% powDiff_perm = squeeze(nanmean(powDiff_perm,1));
% diff_thUP    = prctile(powDiff_perm,95,2);
% diff_thLO    = prctile(powDiff_perm, 5,2);
% powDiff      = nanmean(powDiff,1);
% 
% % visualize
% figure; hold on;
% plot(powDiff, 'linew', 3)
% xlim([0 100])
% plot(diff_thUP, 'r--', 'linew', 1)
% plot(diff_thLO, 'r--', 'linew', 1)
% xlabel('Frequency [HZ]')
% ylabel({'Power Difference Between'; 'Indexed and non-Indexed Trials'})
% mAx = gca;
% mAx.FontWeight = 'bold';
% mAx.FontSize = 18;
% 
% % hz = linspace(0, data.srate, data.pnts);
% % for sub = 1:size(indLFP,1)
% %     subplot(4,2,sub)
% %     plot(hz, nonLFP(sub,:));
% %     hold on
% %     plot(hz, indLFP(sub,:));
% %     xlim([0 50])
% % end