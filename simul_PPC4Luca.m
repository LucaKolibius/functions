% You'll need to add these paths on RDS2 (hanslmas-02); which in my case is
% mapped onto the V: drive
addpath('V:\hanslmas-02\github\Osci_TGMs');
addpath('V:\hanslmas-02\fred4simon\associativeLearningSingleUnit\spectral');
% addpath('V:\hanslmas-02\Common\fieldtrip-20200310\fieldtrip-20200310');

cfg.numchan = 16;
cfg.SR=1000; % Sampling rate
cfg.channel=1;% Number of channels
cfg.time=[-1:0.05:2]; % Time period to be simulated
cfg.trls=100; % Number of trials
cfg.freqs=70; % Frequency of oscillation (will be a sine wave)
cfg.noiseexp=0.6; % Steepness of 1/f 
cfg.amp=0.4; % Amplitude of Oscillation
cfg.jitter = 5;
dataSH=sh_simulEEG_Spks(cfg);

cfg=[];
cfg.channel=1;
EEGdata=ft_selectdata(cfg,dataSH);

cfgtf=[];
cfgtf.method='wavelet';
cfgtf.width=8;
cfgtf.toi=0:0.001:1;
cfgtf.foi=[10:2:150];
cfgtf.output='fourier';
fspec=ft_freqanalysis(cfgtf,EEGdata);
fspec.powspctrm=abs(fspec.fourierspctrm);
powtf=squeeze(nanmean(fspec.powspctrm,1));
ps=squeeze(nanmean(powtf,2));
[tmpf,tmpl]=findpeaks(ps,'SortStr','descend');% find dominant oscillation 

% Extract complex values for PPC calculation at spike times
phi=squeeze(fspec.fourierspctrm(:,1,:,:));% extract phase
ntrls=length(dataSH.trial);
toi(1)=find(dataSH.time{1,1}==cfgtf.toi(1));
toi(2)=find(dataSH.time{1,1}==cfgtf.toi(end));
for nt=1:ntrls
    spkts(nt,:)=dataSH.trial{1,nt}(2,toi(1):toi(2));
end

[ tmp ] = phi;
[ ts ]  = spkts;
ix = find(ts);
tmp1 = permute(phi,[1 3 2]);
tmp2 = reshape(tmp1,[size(tmp1,1)*size(tmp1,2) size(tmp1,3)]);
[ PPC ] = computeSPK2LFPcoupling_S( tmp2(ix,:), 'ppc');
                                                    

%% Plot Powerspectrum of simulated Signal, just to check it looks EEG'ish;
% You can comment this whole section out once you are satisfied with the
% EEG signal
figure;
subplot(1,10,1:6);
pcolor(fspec.time,fspec.freq,squeeze(powtf));
shading interp
ylim([10 150]);
hold on
subplot(1,10,7:8);
plot(ps,fspec.freq);
ylim([10 150]);
hold off
subplot(1,10,9:10);
plot(PPC,fspec.freq);
ylim([10 150]);
hold off

