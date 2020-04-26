function gedSFC_method1(data, idxSpks)
freq    = 3;
chans   = 8;
srate   = 1000;
% nwin    = ceil(srate/freq/8); % window size is 1/4 cycle (1/8 of either side)
nwin = 125;
npnts   = 3000;
trls    = floor(size(data.lfp_noArt,2)/npnts);
covSTrl = zeros(chans);
covS    = zeros(chans);
bbcov   = zeros(chans);
numspks = 0;

for trl = 1:trls
    
%     broadband covariance
    trlLFP = data.lfp_noArt(:,npnts*(trl-1)+1:npnts*trl);
    trlLFP = bsxfun(@minus,trlLFP,mean(trlLFP,2)); % demean broadband data
    
    if any(isnan(trlLFP(:))) % skip if there are artefacts in this trial
        continue
    end
    bbcov  = bbcov + ((trlLFP*trlLFP')/npnts/trls); % compute broadband covariance matrix
    
%     segment spikestimes
    spikelocs = round(idxSpks(idxSpks >= npnts*(trl-1)+1 & idxSpks <= npnts*trl) - npnts*(trl-1));
    spikelocs(spikelocs<=nwin | spikelocs+nwin > npnts) = [];
    
    if isempty(spikelocs) % no spikes within limits of this trial
        continue
    end
    
    numspks = numspks + length(spikelocs);
    % add 4hz noise as a sanity check if the method picks it up
    f = 4;
    t = 0:1/srate:1/f;
    mysine = sin(2*pi*f*t)*20;
    mysine = repmat(mysine,[8,1]);
    trlLFP(:,spikelocs-125:spikelocs+125) = trlLFP(:,spikelocs-125:spikelocs+125) + mysine;

%     spike-locked covariance
    for ti=1:length(spikelocs)
        tempdat = data.lfp_noArt(:,spikelocs(ti)-nwin:spikelocs(ti)+nwin); % take the EEG activity from all channels at timestamps surrounding the first trough
        tempdat = bsxfun(@minus,tempdat,mean(tempdat,2)); % demean data snippet
        covSTrl   = covSTrl + (tempdat*tempdat')/nwin; % calculate covariance of trough locked snippet and add it to existing covariance matrix
    end
    covS = covS + (covSTrl./ti/trls); % average by number of added troughs
    
end


[evecs,evals] = eig(covS,bbcov);
[~,compidx]   = sort(diag(evals)); % max component


for trl = 1 : trls
    trlLFP = data.lfp_noArt(:,npnts*(trl-1)+1:npnts*trl);
    trlLFP = trlLFP'*evecs(:,compidx(end));
    
    trlPow(trl,:) = abs(fft(trlLFP/npnts));
end

%% Visualisation
% _______________

% vector of frequencies
hz = linspace(0,srate,npnts);

figure
% subplot(331)
plot(hz, nanmean(trlPow,1));
set(gca,'xlim',[0 90],'xtick',10:20:90), axis square
xlabel('Frequency (Hz)'), ylabel('Power (a.u.)')
title(sprintf('%d# Spikes', numspks));
end


%% for continuous data:
% freq    = 4;
% chans   = 8;
% srate   = 1000;
% nwin    = ceil(srate/freq/8); % window size is 1/4 cycle (1/8 of either side)
% data    = data.lfp_raw;
% npnts   = size(data,2);
% trls    = 1;
% covSTrl = zeros(chans);
% covS    = zeros(chans);
% bbcov   = zeros(chans);
% 
% for trl = 1:trls
%     
%     % broadband covariance
%     trlLFP = data(:,npnts*(trl-1)+1:npnts*trl);
%     trlLFP = bsxfun(@minus,trlLFP,mean(trlLFP,2)); % demean broadband data
%     
%     if any(isnan(trlLFP(:))) % skip if there are artefacts in this trial
%         continue
%     end
%     bbcov  = bbcov + ((trlLFP*trlLFP')/npnts/trls); % compute broadband covariance matrix
%     
%     % segment spikestimes
%     spikelocs = round(idxSpks(idxSpks >= npnts*(trl-1)+1 & idxSpks <= npnts*trl) - npnts*(trl-1));
%     spikelocs(spikelocs<=nwin | spikelocs+nwin > npnts) = [];
%     
%     if isempty(spikelocs) % no spikes within limits of this trial
%         continue
%     end
%     
%     % spike-locked covariance
%     for ti=1:length(spikelocs)
%         tempdat = data(:,spikelocs(ti)-nwin:spikelocs(ti)+nwin); % take the EEG activity from all channels at timestamps surrounding the first trough
%         tempdat = bsxfun(@minus,tempdat,mean(tempdat,2)); % demean data snippet
%         covSTrl   = covSTrl + (tempdat*tempdat')/nwin; % calculate covariance of trough locked snippet and add it to existing covariance matrix
%     end
%     covS = covS + (covSTrl./ti/trls); % average by number of added troughs
%     
% end
% 
% 
% [evecs,evals] = eig(covS,bbcov);
% [~,compidx]   = sort(diag(evals)); % max component
% 
% 
% for trl = 1 : trls
%     trlLFP = data(:,npnts*(trl-1)+1:npnts*trl);
%     trlLFP = trlLFP'*evecs(:,compidx(end));
%     
%     trlPow(trl,:) = abs(fft(trlLFP/npnts));
% end

