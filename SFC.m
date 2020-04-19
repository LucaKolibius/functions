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

it = 2;
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

%%
% get spkDat
% get all spks from that session
load(['X:\Luca\data\allSbj\allSU_', bidsID,'_', sesh,'.mat'], 'allSU')

IUspks = allSU{su,3}/32; % get from 32khz to 1khz

% split spikes into spikes that occur in index trials and those that do not
idxTW = data.tw(data.idxTrl,:); % time windows of index trials
for ia = 1:size(idxTW,2)
    idx = IUspks>=idxTW(ia,1) & IUspks<=idxTW(ia,2);
    spksInd   = round(IUspks ( idx )); % spiketimes (1 kHz) during     indexed trials
    spksNoInd = round(IUspks (~idx )); % spiketimes (1 kHz) during non-indexed trials
end

% lfp dat is in data.lfp_noArt
npad   = 40; % even only, please!
npad2  = npad/2;
npnts  = size(data.lfp_noArt,2);
nchans = 8;
srate  = 1000;

% produce augmented data (delay embedding)
padorder = [ npnts-floor(npad2):npnts 1:floor(npad2)-1 ];

delEmb = zeros(nchans*npad,npnts);
for deli = 1:npad
    delEmb( (1:nchans)+(deli-1)*nchans,:) = detrend(data.lfp_raw(:,[padorder(deli):end 1:padorder(deli)-1])')';
end

% sphere data by scaling the eigenvectors
[evecsO,evalsO] = eig( (delEmb*delEmb')/size(delEmb,2) );
spheredata = (delEmb' * evecsO * sqrt(inv(evalsO)) )';
spheredata = reshape(spheredata,[npad*nchans npnts ]);


% sum covariances around spikes, then divide by N
% spikelocs is an index (of 500 spikes)
spcov = zeros(size(delEmb,1));
% spikelocs = spksInd;
spikelocs = round(IUspks);
for si=1:length(spikelocs)
    if spikelocs(si)-npad < 0 || spikelocs(si)+npad > size(spheredata,2)
        continue
    end
    tmpdat = spheredata(:,spikelocs(si)-npad:spikelocs(si)+npad);
    spcov  = spcov + tmpdat*tmpdat'/size(tmpdat,2);
end
spcov = spcov/si;


% eigendecomposition of sphered matrix
[evecsF,evalsF] = eig( spcov );

% compute weights and map
% you here rotate the initial
jdw    = evecsO * sqrt(pinv(evalsO)) * evecsF;
jdmaps = pinv(jdw)';


% simple spike-triggered average
sta = zeros(nchans,npad+1);
for si=1:length(spikelocs)
    if spikelocs(si)-npad2 < 0 || spikelocs(si)+npad2 > size(data.lfp_raw,2)
        continue
    end
    sta = sta+data.lfp_raw(:,spikelocs(si)-npad2:spikelocs(si)+npad2);
end

%% plotting

figure(1), clf
tv = 1000*(-npad2:npad2)/srate;
addWire = 6;


% forward model of spatiotemporal component
subplot(221)
rmap = reshape(jdmaps(:,end)',nchans,npad);
contourf(npad2+tv(1:end-1),1:nchans,rmap,40,'linecolor','none')

hold on, plot(npad2+tv(1:end-1),4+1*zscore(rmap(wireU,:)),'k','linew',2)
hold on, plot(npad2+tv(1:end-1),4+1*zscore(rmap(addWire,:)),'r','linew',2)
hold on, plot([20 20],get(gca,'ylim'),'k--')


title('Spatiotemporal component'), %axis square
xlabel('Filter time (ms)'), ylabel('"Cortical depth"')

% Spike-triggered average
subplot(222)
imagesc(sta)
contourf(tv,1:nchans,sta,40,'linecolor','none')
hold on, plot([0 0],get(gca,'ylim'),'k--')
% staDat = 4*sta(6,:)./max(sta(6,:));
% hold on; plot(staDat, 'linew', 2, 'color', 'k');
title('Spike-triggered average'), %axis square
xlabel('Peri-spike time (ms)'), ylabel('"Cortical depth"')

% power spectrum of component and STA
subplot(2,4,[5 6])
plot(linspace(0,srate,200),abs(fft(zscore(rmap(wireU,:))/npad,200)).^2,'ks-','linew',2,'markersize',8,'markerfacecolor','w'), hold on
plot(linspace(0,srate,200),abs(fft(zscore(sta(wireU,:)) /npad,200)).^2,'bo-','linew',2,'markersize',8,'markerfacecolor','w')
set(gca,'xlim',[0 100])
xlabel('Frequencies (Hz)'), ylabel('Power (a.u.)')
title(sprintf('MW %d (IU wire)', wireU))

subplot(2,4,[7 8])
plot(linspace(0,srate,200),abs(fft(zscore(rmap(addWire,:))/npad,200)).^2,'ks-','linew',2,'markersize',8,'markerfacecolor','w'), hold on
plot(linspace(0,srate,200),abs(fft(zscore(sta(addWire,:)) /npad,200)).^2,'bo-','linew',2,'markersize',8,'markerfacecolor','w')
set(gca,'xlim',[0 100])
xlabel('Frequencies (Hz)'), ylabel('Power (a.u.)')
legend({'STF';'STA'})
title(sprintf('MW %d (suppl. wire)', addWire))

print(gcf,['X:\Luca\SFC\SFC_', bidsID, '_', sesh, '_SU', num2str(su),'.jpg'],'-djpeg','-r720');
