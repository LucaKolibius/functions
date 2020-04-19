% This function simulates a simple EEG with one ongoing oscillation and
% spikes that are coupled to that oscillation
% overlaid by a 1/f type noise; Output is in a fieldtrip type structure


function data=sh_simulEEG_Spks(cfg)
numchan = cfg.numchan;
Fs=cfg.SR; % sample rate
dt=1/Fs;
x=cfg.time(1):dt:cfg.time(end); % timevector
elength=numel(x); % experiment length
ntrls=cfg.trls; % number of trials
f=cfg.freqs; % frequency of interest
jitter=cfg.jitter;% defines the maximal jitter of spike times in ms; Spikes will be uniformly distributed in that time window
exp=cfg.noiseexp;
amp=cfg.amp;
dumspks=zeros(1,elength);
data.label={'EEG';'Spks'};
for n=1:cfg.trls
    % Simulate a simple signal
    % use different values for exp to create different 1/f shapes
    for nn = 1:numchan
    noise(nn,:)=pinknoise(length(x),exp);
    end 
    
    % jitter phase accross trials
    phi=rand(1,1)*2*pi; % randomly jitter the sine wave for each trial
    tmpsig = sin(2*pi*x*f+phi); % generate sine wave 
    [~,locs]=findpeaks(tmpsig); % find peak in the sine wave
    jitt=(rand(1,numel(locs))-0.5).*jitter; % put the spikes around the peaks (locs) with a jitter
    spkidx=round(locs+jitt);
    
    % a quick fix
    if spkidx(1) == 0
        spkidx(1) = [];
    end
    tmpsig = repmat(tmpsig, [numchan, 1]);
    
    trunc=find(spkidx>elength);
    if ~isempty(trunc)
        spkidx(trunc)=[];
    end
    tmpspks=dumspks;
    tmpspks(1,spkidx)=1;
    y(1,:) = tmpsig(1,:).*amp+noise(1,:);
    y(2,:) = tmpspks;
    data.trial{1,n}=y; 
    data.time{1,n}=x;
    data.sampleinfo(n,1)=(n-1)*length(x)+1;
    data.sampleinfo(n,2)=data.sampleinfo(n,1)+length(x)-1;
    data.trials{1,n} = tmpsig.*amp+noise;
    data.trials{2,n} = tmpspks;
end
end
