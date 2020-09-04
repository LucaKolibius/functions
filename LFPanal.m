% various settings
sets         = [];
sets.server  = 0;    % 0 for local machine, 1 for server
sets.pddng   = 500;  % padding before and after data snippet for SFC to increase frequency resolution
sets.maxSpks = 5000; % maximal number of spikes that is considered for SFC
sets.srate   = 1000; % sampling rate
sets.nwin    = 50;   % time in samples before and after spikes that is considered for analysis
sets.chan    = 8;    % number of EEG channels
sets.nperm   = 1000; % number of permutations

switch sets.server
    case 0 % local machine
        addpath('X:\Luca\toolboxes\gedCFC_tutorial');
        IUdir  = dir('X:\Luca\indexSUvisu\varWindow_ENC_-1-113_RET_-1113_DP_p99_th1\SU_*.mat');
        LFPdir = dir('X:\George\Analysis\Data Continuous\sub-*_micro_RAW_DS-1000_SPK-INT.mat');
        artDir = dir('X:\Luca\engram allocation lfp\lfpArtefacts\lfp_artefacts\LFPartefacts_*.mat');
    case 1 % server
        addpath('/media/ldk898/rds-share/Luca/server/functionen');
        IUdir  = dir('/media/ldk898/rds-share/Luca/server/varWindow_ENC_-1-113_RET_-1113_DP_p99_th1/SU_*.mat');
        LFPdir = dir('/media/ldk898/rds-share/Luca/server/George/sub-*_micro_RAW_DS-1000_SPK-INT.mat');
        artDir = dir('/media/ldk898/rds-share/Luca/server/lfp_artefacts/LFPartefacts_*.mat');
end

% start with each SU
for it = 1:size(IUdir,1)
    clearvars -except IUpath IUdir LFPdir artDir it nonLFP_all indLFP_all sets permLFP numIdx powDiff logThis
    switch sets.server
        case 0
            load([IUdir(it).folder, filesep, IUdir(it).name]) % load IU info
            load([artDir(1).folder, filesep, 'LFPartefacts_', bidsID, '_', sesh, '_enc_', wireName(1:end-1),'.mat'], 'artefactsEnc'); % load artefact timestamps (1x8 cell with each cell corresponding to a wire)
        case 1
            load(['/media/ldk898/rds-share/Luca/server/varWindow_ENC_-1-113_RET_-1113_DP_p99_th1', '/', IUdir(it).name]) % load IU info
            load(['/media/ldk898/rds-share/Luca/server/lfp_artefacts', '/', 'LFPartefacts_', bidsID, '_', sesh, '_enc_', wireName(1:end-1),'.mat'], 'artefactsEnc'); % load artefact timestamps (1x8 cell with each cell corresponding to a wire)
    end
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
    switch sets.server
        case 0
            cd(['X:\Luca\data\', mSubject,filesep,  mSession,filesep]);
            abc  = dir('2*'); cd(abc.name);
            p2d  = [cd, filesep];
            trls = 1; % can be ignored here; for when you want to have a specific animal cue
            [~, hitsIdx, ~, ~, ~, ~, encTrigger, ~, ~] = loadLogs(p2d, trls);
        case 1
            load(['/media/ldk898/rds-share/Luca/server/logs/', num2str(it)])
    end
    
    % encoding trial
    encTrigger = encTrigger(hitsIdx,1);  % only hits, only cue-locked
    encTrigger = round(encTrigger*1000); % seconds to samples
    
    % load wirename
    switch sets.server
        case 0
            loadFile = dir(['X:\Luca\data\allSbj\allSU_', bidsID,'_',  mSession '.mat']);
            load([loadFile.folder, filesep, loadFile.name])
        case 1
            loadFile = dir(['/media/ldk898/rds-share/Luca/server/data/allSbj/allSU_', bidsID,'_',  mSession '.mat']);
            load(['/media/ldk898/rds-share/Luca/server/data/allSbj', '/', loadFile.name])
    end
    
    idxWire   = allSU{su,1};
    idxNum    = str2double(idxWire(end));
    idxBundle = idxWire(1:end-1);
    disp(idxWire);
    
    loadLFP = cellfun(@(x)regexp(x,[bidsID, '_', mSession]), {LFPdir.name}, 'UniformOutput', false); % find the LFP date for the correct patient and session ...
    loadLFP = find(~cellfun(@isempty, loadLFP));
    
    switch sets.server
        case 0
            load([LFPdir(loadLFP).folder, filesep, LFPdir(loadLFP).name]) % ... and load it in
        case 1
            load(['/media/ldk898/rds-share/Luca/server/George', '/', LFPdir(loadLFP).name]) % ... and load it in
    end
    
    idx = cellfun(@(x)regexp(x,idxBundle), data_micro.label, 'UniformOutput', false); % find the position of the MW of my bundle of interest
    idx = cellfun(@(x)find(x==1), idx, 'UniformOutput', false);
    idx = ~cellfun(@isempty, idx);
    idx = find(idx == 1);
    
    data = [];
    data.lfp_raw = data_micro.trial{1}(idx,:);
    data.label   = data_micro.label(idx);
    data.idxTrl  = sensTrlsFFPP;
%     data.tw      = [encTrigger-1000 encTrigger];
    data.tw      = [encTrigger-2000 encTrigger+1000];
    
    % reject artefacts
    for artMW  = 1 : size(data.lfp_raw,1) % for each MW
        artVec = []; % generate a vector with all the timestamps at which we have an artefact
        for t = 1:size(artefactsEnc{artMW},1) % create a vector with timepoints of all artefact sample points
            artVec = [artVec artefactsEnc{artMW}(t,1) : artefactsEnc{artMW}(t,2)];
        end
        tempLFP                 = data.lfp_raw(artMW,:);
        tempLFP(artVec)         = NaN; % make these timepoints NaNs
        data.lfp_noArt(artMW,:) = tempLFP; % save this LFP as artefact rejected LFP (lfp_noArt)
    end
    
    % next we look at trials and if they contain more than a predefined number
    % of nans. If they do, we chuck the whole trial (we don't do that
    % anymore)
    for trl = 1:size(encTrigger,1)
        segmDat                = data.lfp_noArt(:,data.tw(trl, 1):data.tw(trl, 2));
        segmDatShort           = data.lfp_noArt(:,encTrigger(trl)-550:encTrigger(trl)+550); % has to be 1101 long;
        %         delTrl            = sum(isnan(segmDat),2) > 500; % more than half of the segment is an artefact
        %         segmDat(delTrl,:) = NaN; % chuck the whole segment
        data.segmDat{trl}      = segmDat;
        data.segmDatShort{trl} = segmDatShort;
    end
    
%     data.pnts  = size(data.segmDat{1},2);
%     data.srate = sets.srate;
%     data.chan  = size(data.label,1);
    
    % Get the spikes
    % if the number of spikes exceed sets.maxSpks, draw a number of random
    % spikes specified by sets.maxSpks
    idxSpks = allSU{su,3}/32;
    if length(idxSpks) > sets.maxSpks
        randDraw = randperm(length(idxSpks));
        randDraw = randDraw(1:sets.maxSpks);
        randDraw = sort(randDraw, 'ascend');
        idxSpks = idxSpks(randDraw);
    end
    
    %% load results if they are already there
    SFCout = dir('X:\Luca\engram allocation lfp\SFC\method1\SFC_*');
    
    % if the SFC already was calculated, load it
    if exist([SFCout(1).folder, filesep, 'SFC_', bidsID, '_', sesh, '_su', num2str(su),'.mat'], 'file') 
        % load eigenvector of interest "evec" and check if it signfiicant coupling
        load([SFCout(1).folder, filesep, 'SFC_', bidsID, '_', sesh, '_su', num2str(su),'.mat'], 'evec', 'signSFC', 'evecs')   
        
        % check to which frequencies the IU couples
        figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        hzPow = linspace(0,sets.srate, (sets.nwin*2)+1 + 2*sets.pddng);
        
        load([SFCout(1).folder, filesep, 'SFC_', bidsID, '_', sesh, '_su', num2str(su),'.mat'], 'trlPow', 'trlPowPerm')
        tmp  = evec == evecs; comp = find(tmp(1,:) == 1);
        
        trlPow = squeeze(trlPow(comp,:,:));
        trlPow = nanmean(trlPow,1);
        
        % plot thresholds
        plot(hzPow, prctile(trlPowPerm, 95,1), 'color', 'g'); hold on
        plot(hzPow, prctile(trlPowPerm,  5,1), 'color', 'g');
        
        TSsig = trlPow;
        TSnon = trlPow;
        
        TSsig(trlPow < prctile(trlPowPerm,95,1))  = nan;
        TSnon(trlPow >= prctile(trlPowPerm,95,1)) = nan;
        
        plot(hzPow, TSsig, 'r', 'linew', 2);
        plot(hzPow, TSnon, 'color', [0 0 0], 'linew', 2);
        xlim([1 250])
        
        sigSFCfreq = hzPow(trlPow>=prctile(trlPowPerm,95,1)); % these are the frequency to which the IU significantly locks
    else
        break % they should all be there already
        filtf   = 0; % highpass filter at 10hz; set to 0 for no filtering
        [evec, signSFC] = gedSFC_method1(data, idxSpks, filtf, bidsID, sesh, su, sets);
        
        if signSFC == 0
            filtf   = 10; % highpass filter at 10hz; set to 0 for no filtering
            [evecf, signSFC] = gedSFC_method1(data, idxSpks, filtf, bidsID, sesh, su, sets); % redo SFC
            
            % update evec if the filtered approach yielded significant SFC
            if signSFC == 1
                evec = evecf;
            end
        end
    end
    
    %     normal, unweighted average
    %     evec = ones(8,1)*0.125; % normal average
    
    %     only the wire on which the SU is found
    %     evec = zeros(8,1);
    %     evec(idxNum) = 1;
    
    %   only the wire with the highest scalar in evec
    [~,b]= max(evec);
    evec = zeros(8,1);
    evec(b) = 1;
    
    logThis(1,it) = idxNum; % SU found on wire idxNum
    logThis(2,it) = b;      % max SFC
    
    data.segmComp      = cellfun(@(x) x'*evec, data.segmDat, 'un', 0);
    data.segmCompShort = cellfun(@(x) x'*evec, data.segmDatShort, 'un', 0);
    
    % for time frequency
    indLFP  =  horzcat(data.segmComp{ data.idxTrl});
    nonLFP  =  horzcat(data.segmComp{~data.idxTrl});
    numIdx  =  sum(data.idxTrl);
    
    % for power we take a shorter segment to match the SFC
    indLFPshort = horzcat(data.segmCompShort { data.idxTrl});
    nonLFPshort = horzcat(data.segmCompShort {~data.idxTrl});
    
    switch sets.server
        case 0
            %             savepath = ['X:\Luca\engram allocation lfp\SFC\method1\', num2str(filtf), filesep, num2str(sets.nwin), filesep];
            savepath = 'X:\Luca\engram allocation lfp\SFC\method1\';
        case 1
            %             savepath = ['/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/fileshare/engram allocation lfp/', num2str(filtf), '/', num2str(sets.nwin), '/'];
            savepath = ['/media/ldk898/rds-share/Luca/server/engram allocation lfp/', num2str(filtf), '/', num2str(sets.nwin), '/'];
    end
    save([savepath,  'LFP_', bidsID, '_', sesh, '_su', num2str(su)], 'numIdx', 'indLFP', 'nonLFP', 'sigSFCfreq', 'indLFPshort', 'nonLFPshort')
    
end

%% powDiff
clear
clc
LFPidx   = dir('X:\Luca\engram allocation lfp\SFC\method1\LFP*');
numIdxA  = [];
datpnts  = 1101;
for iu = 1 : size(LFPidx,1)
    load([LFPidx(iu).folder, filesep, LFPidx(iu).name], 'indLFPshort', 'nonLFPshort', 'numIdx', 'sigSFCfreq')
    indLFP = indLFPshort;
    nonLFP = nonLFPshort;
    
    %     % find power spectrum difference per MW and then average
    powInd           = abs(fft(indLFP) / datpnts);
    powNon           = abs(fft(nonLFP) / datpnts);
    powDiff(iu,:)    = nanmean(powInd,2) - nanmean(powNon, 2);
    
    % create a variable that  has all the indLFP and nonLFP powerspectra in one double, separate for each MW
    permLFP{iu}      = [powNon powInd];
    
    numIdxA          = [numIdxA numIdx];
    sigSFCfreqA{iu}  = sigSFCfreq;
end
numIdx  = numIdxA;

clearvars -except numIdx permLFP powDiff LFPidx sigSFCfreqA

nperm = 10000;
for perm = 1:nperm % repeat 10.000 times
%     disp(perm)
    for iu = 1:size(LFPidx,1) % cycle through each IU
        
        numIdxTrl               = numIdx(iu); % number of trials this unit indexes
        permDx                  = randperm(size(permLFP{iu},2));
        powDiff_perm(iu,:,perm) = nanmean(permLFP{iu}(:, permDx(1:end-numIdxTrl)),2) - nanmean(permLFP{iu}(:, permDx(end-numIdxTrl+1:end)), 2);
        
    end
    
end


for iu = 1 : size(LFPidx,1)
    
powDiff_perm_thUP5 = prctile(powDiff_perm(iu,:,:), 99, 3);
powDiff_perm_thLO5 = prctile(powDiff_perm(iu,:,:),  1, 3);

% show significant power differences
figure('units','normalized','outerposition',[0 0 1 1])
hold on
hz = linspace(0,1000,1101);
plot(hz, powDiff(iu,:), 'linew', 2, 'color', [0 0 0]);
plot(hz, powDiff_perm_thUP5, 'color', 'r');
plot(hz, powDiff_perm_thLO5, 'color', 'r');

% mark significant spike field coupling as a red bar
sigSFCfreq  = sigSFCfreqA{iu};
tmp = ~ismember(hz, sigSFCfreq);
hz(tmp) = nan;

posi = get(gca, 'Ylim');
plot([hz], ones(1,1101)*posi(2), 'linew', 5, 'color', 'r')
xlim([1 200])

aa (1,iu) = mean(powDiff(iu,~tmp));
aa (2,iu) = mean(powDiff_perm_thUP5(~tmp));
drawnow
saveas(gcf, ['X:\Luca\engram allocation lfp\results\powerDiff_evecWire\powerDiff_evecWire__', num2str(iu), '.png'])

end

%% now look at the TF
% make variable to fieldtrip struct
LFPidx      = dir('X:\Luca\engram allocation lfp\SFC\method1\LFP*');
dat.fsample = 1000;
dat.npnts   = 3001;
dat.label   = {'chan1'};
close all
for iu = 1 : size(LFPidx,1)
    load([LFPidx(iu).folder, filesep, LFPidx(iu).name], 'indLFP', 'nonLFP', 'numIdx', 'sigSFCfreq')
    
    datNon = dat;
    for trl = 1 : size(nonLFP,2)
        datNon.trial{trl} = nonLFP(:,trl)';
        datNon.time {trl} = 1/dat.fsample : 1/dat.fsample : dat.npnts/dat.fsample;
    end
    
    datInd = dat;
    for trl = 1 : size(indLFP,2)
        datInd.trial{trl} = indLFP(:,trl)';
        datInd.time {trl} = 1/dat.fsample : 1/dat.fsample : dat.npnts/dat.fsample;
    end
    
    cfgtf           = [];
    cfgtf.method    = 'wavelet';
    cfgtf.width     = 5;
    cfgtf.toi       = 'all';
    cfgtf.foi       = [1:1:200];
    cfgtf.output    = 'fourier';
    fspec           = ft_freqanalysis(cfgtf,datInd);
    tfInd           = fspec.fourierspctrm;
    
    fspec           = ft_freqanalysis(cfgtf,datNon);
    tfNon           = fspec.fourierspctrm;
    
    % time, frequency, trial, chanel
    % trials, channel, freq, time
    figure('units','normalized','outerposition',[0 0 1 1])
    lims = [];
    for trl = 1 : size(tfInd,1)
        subplot(3, 3, trl)
        plotThis = squeeze(tfInd(trl,:,:,:));
        imagesc(abs(plotThis));
        lims(trl,:) = caxis;
        title('Indexed')
        xlim([1000 2000])
        xticks(1000:200:2000);
        xticklabels([0:200:1000])
        hand = gca;
        hand.YDir = 'normal';
%         colorbar
    end
    
    for trl = 1 : size(tfInd,1)
        subplot(3, 3, trl)
        caxis(gca, [min(lims(:,1)) max(lims(:,2))/10 ]);
% caxis(gca, [0 50 ]);
    end
        
        
    for trl = 1 : 9-size(tfInd,1)
        subplot(3, 3, size(tfInd,1) + trl)
        plotThis = squeeze(tfNon(trl,:,:,:));
        imagesc(abs(plotThis));
                caxis(gca, [min(lims(:,1)) max(lims(:,2))/7 ]);
%         caxis(gca, [0 50 ]);
        title('non-indexed');
        xlim([1000 2000])
        xticks(1000:200:2000);
        xticklabels([0:200:1000])
        hand = gca;
        hand.YDir = 'normal';
        %         colorbar
    end
    drawnow
    saveas(gcf, ['X:\Luca\engram allocation lfp\results\powerDiff_evecWire\tf_power\tf_pow_evecWire__', num2str(iu), '.png'])

    % angle
    figure('units','normalized','outerposition',[0 0 1 1])
    test = fspec.fourierspctrm(trl,1,:,:);
    subplot(1,20,1:19)
test = squeeze(test);
test = angle(test);
imagesc(test)
xlim([1000 2000])
xticks(1000:200:2000);
xticklabels([0:200:1000])
xlabel('Time in ms')
ylabel('Frequency in hertz')
hand = gca;
hand.FontWeight = 'bold';
hand.FontSize = 20;

subplot(1,20,20)
pow = abs(squeeze(fspec.fourierspctrm(1,1,:,:)));
tem = nanmean(pow(:,1001:2000), 2);
plot(tem, 'linew', 2, 'color', 'k');
view([90 90])
xticks('')
yticks('')
drawnow

saveas(gcf, ['X:\Luca\engram allocation lfp\results\powerDiff_evecWire\tf_phase\tf_phs_evecWire__', num2str(iu), '.png'])

end





%% old stuff
%% OLD SFC
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

%
% this part is the old poweranalysis of the 1s before the cue
%     numtrl = size(data.idxTrl,2);
%     nonLFP = []; indLFP = [];
%
%     nonCounter = 0; % start at trial 1
%     indCounter = 0;
%
%     for trl = 1:numtrl
%         if data.idxTrl(trl) == 0 % not indexed
%             nonCounter = nonCounter + 1;
%             nonLFP(nonCounter,:) = abs(fft(data.segmComp{trl}) / data.pnts);
%         elseif data.idxTrl(trl) == 1 % indexed trial
%             indCounter = indCounter + 1;
%             indLFP(indCounter,:) = abs(fft(data.segmComp{trl}) / data.pnts);
%         end
%     end
%
%     %% find empirical difference
%     % find power spectrum difference per MW and then average
%     powDiff(it,:) = nanmean(nonLFP,1) - nanmean(indLFP, 1);
%
%     % create a variable that  has all the indLFP and nonLFP powerspectra in one
%     % double, separate for each MW
%
%     permLFP{it} = [nonLFP; indLFP];
