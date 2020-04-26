%% multivariate cross-frequency coupling, based on generalized eigendecomposition (gedCFC)
% Method 5: spike-field coherence adapted from Methods 3 and 4 (sphered delay-embedded matrices)

% You will need the following files in the current directory or Matlab path:
%   - filterFGx.m

% mikexcohen@gmail.com

% Adapted for multiple triasl by LDK 23.04.2020
% this script has been used to analyze simulated data from an adapted
% simulation script from Simon

%% the important stuff...
% we need to get this from simons simulation
npad = 40;
npad2 = npad/2;
nchans = 16;
srate = 1000;
npnts  = 3001;
trls = 100;

% % produce augmented data
padorder = [ npnts-floor(npad2):npnts 1:floor(npad2)-1 ];
delEmb = zeros(nchans*npad,npnts);
delEmbCov = zeros(nchans*npad); 
spcovTrl = zeros(nchans*npad);
for trl = 1:trls % loop over trials
    data = dataSH.trials{1,trl}(:,:); %% trial data

for deli = 1:npad % for this trial, compute delay embedding
    delEmb( (1:nchans)+(deli-1)*nchans,:) = detrend(data(:,[padorder(deli):end 1:padorder(deli)-1])')'; % 9980 until end of data + 1 until 9979
end
trialCov = delEmb*delEmb';

delEmbTrl{trl} = delEmb;
delEmbCov = delEmbCov + trialCov / trls;
end

% sphere data
[evecsO,evalsO] = eig( delEmbCov / size(delEmbCov,2) );

% lets try to sphere each trial
for trl = 1:trls
    tmp = (delEmbTrl{trl}' * evecsO * sqrt(inv(evalsO)) )';
    spheredata{trl} = reshape(tmp,[npad*nchans npnts ]);
end

% sum covariances around spikes, then divide by N
for trl = 1:trls
    spcov = zeros(size(delEmbTrl,1));
    
    % get spiketimes
    % spikes
    spikelocs = dataSH.trials{2,trl};
    spikelocs = find(spikelocs == 1); % spike position
    % quick fix in case a spike is outside the time vector
    spikelocs(spikelocs-npad <= 0) = [];
    spikelocs(spikelocs+npad > npnts) = [];
    
    for si=1:length(spikelocs)
        tmpdat = spheredata{trl}(:,spikelocs(si)-npad:spikelocs(si)+npad);
        spcov  = spcov + tmpdat*tmpdat'/size(tmpdat,2);
    end
    
    spcov = spcov/si; % normalize by number of spikes
    spcovTrl = spcovTrl + spcov / trls;
end

% eigendecomposition of sphered matrix
[evecsF,evalsF] = eig( spcovTrl );

% compute weights and map
jdw    = evecsO * sqrt(pinv(evalsO)) * evecsF;
jdmaps = pinv(jdw)';

%% plotting

figure(1), clf
tv = 1000*(-npad2:npad2)/srate;

% forward model of spatiotemporal component
subplot(121)
rmap = reshape(jdmaps(:,end)',nchans,npad);
contourf(npad2+tv(1:end-1),1:nchans,rmap,40,'linecolor','none')
hold on, plot(npad2+tv(1:end-1),8+1*zscore(rmap(8,:)),'k','linew',2)
title('Spatiotemporal component'), axis square
xlabel('Filter time (ms)'), ylabel('"Cortical depth"')

% power spectrum of component
subplot(122)
plot(linspace(0,srate,200),abs(fft(zscore(rmap(1,:))/npad,200)).^2,'ks-','linew',2,'markersize',8,'markerfacecolor','w'), hold on
set(gca,'xlim',[0 200])
xlabel('Frequencies (Hz)'), ylabel('Power (a.u.)')
legend('STF')

%%
plot(abs(fft(data(1,:))))
set(gca,'xlim',[0 200])