function [evec, signSFC] = gedSFC_method1(data, idxSpks, filtf, bidsID, sesh, su, sets)
% find empirically optimal time windows for segments
[tw1, tw2] = fndLFPsegs(data, idxSpks, sets.nwin);


npnts   = tw2(1)-tw1(1); % they are all the same length
trls    = size(tw1,1);

chans   = sets.chan;
srate   = sets.srate;
covS    = zeros(chans);
bbcov   = zeros(chans);
numspks = 0;
effTrls = 0;
spkLFP  = [];
nperm   = sets.nperm;
pddng   = sets.pddng; % padding
sniplen = ((sets.nwin*2)+1) + 2*pddng; % snippet length


for trl = 1:trls
    
    %     broadband covariance
    trlLFP = data.lfp_noArt(:,tw1(trl) : tw2(trl));
    
    %% highpass filter data (doing the filtering on the snippits is suboptimal)
    if filtf ~= 0
        trlLFP = eegfilt(trlLFP, srate, filtf, 0);
    end
    
    trlLFPd = bsxfun(@minus,trlLFP,mean(trlLFP,2)); % demean broadband data
    
    % spikestimes in segment
    %     spikelocs = round(idxSpks(idxSpks >= npnts*(trl-1)+1 & idxSpks <= npnts*trl) - npnts*(trl-1)); % old
    spikelocs = round(idxSpks(idxSpks >= tw1(trl) & idxSpks <= tw2(trl)) - tw1(trl));
    spikelocs(spikelocs<=sets.nwin | spikelocs+sets.nwin > npnts) = [];
    
    if any(isnan(trlLFPd(:))) % skip if there are artefacts in this trial
        continue
    end
    
    if isempty(spikelocs) % no spikes within limits of this trial
        continue
    end
    
    effTrls = effTrls + 1; % effective number of trials (trials without an artifact
    numspks = numspks + length(spikelocs); % number of spikes used for SFC
    disp(numspks)
    % save spikelocs and LFP for later permutation
    allTrlLFP {effTrls} = trlLFP;
    allTrlSPK {effTrls} = spikelocs;
    
    % calculate broadband covariance matrix
    bbcov  = bbcov + ((trlLFPd*trlLFPd')/npnts+1);
    
    % %     %     % add 20hz noise as a sanity check if the method picks it up
    %         f = 80;
    %         t = 0:0.001:0.02; % sr is 0.001 and duration is 0.1
    %         mysine = sin(2*pi*f*t)*200;
    %         mysine = repmat(mysine,[8,1]);
    %         trlLFP(:,spikelocs-sets.nwin:spikelocs+sets.nwin) = trlLFP(:,spikelocs-sets.nwin:spikelocs+sets.nwin) + mysine;
    
    % I'll go here for the approach to extract the LFP surrounding each spike. There might be an overlap between those time windows (for bursts) and selecting a window in which there are multiple bursts might be a valid alternative. But as I am interested in single spikes here and not bursts, I will consider the LFP surrounding each spike
    for spk = 1:length(spikelocs)
        spkLFP = [spkLFP; {trlLFP(:,spikelocs(spk)-sets.nwin:spikelocs(spk)+sets.nwin)}];
    end
    
    %     spike-locked covariance
    covSTrl = zeros(chans);
    for ti=1:length(spikelocs)
        tempdat = trlLFP(:,spikelocs(ti)-sets.nwin:spikelocs(ti)+sets.nwin); % take the EEG activity from all channels at timestamps surrounding the first spike
        tempdat = bsxfun(@minus,tempdat,mean(tempdat,2)); % demean data snippet
        covSTrl = covSTrl + (tempdat*tempdat')/((sets.nwin*2)+1); % calculate covariance of trough locked snippet and add it to existing covariance matrix
    end
    covS = covS + covSTrl; % average by number of added troughs
    
end
bbcov = bbcov ./ effTrls;
covS  = covS  ./ numspks;

[evecs,evals] = eig(covS,bbcov);
[~,compidx]   = sort(diag(evals)); % max component

trlPow = zeros(chans, numspks, sniplen);
for cmp = chans : -1: 1
    for spk = 1 : size(spkLFP,1)
        spkLFPlin = spkLFP{spk};
        spkLFPlin = spkLFPlin'*evecs(:,compidx(chans));
        spkLFPlin = [zeros(pddng,1); spkLFPlin; zeros(pddng,1)];
        trlPow(cmp, spk,:) = abs(fft(spkLFPlin)/sniplen);
    end
end

% permutation test
trlPowPerm = zeros(nperm, sniplen); % preallocation
hz         = linspace(0, srate, sniplen);
signSFC    = 0;

for cmp = chans : - 1 : 1
    if signSFC == 1
        continue
    end
    tic
    parfor perm = 1:nperm
        covS    = zeros(chans);
        
        % shuffle the spikes
        randex        = randperm(size(allTrlSPK,2));
        allTrlSPKperm = allTrlSPK(randex);
        spkLFPperm    = [];
        
        for permtrl = 1:effTrls
            spikelocs = allTrlSPKperm{permtrl};
            
            % spike-locked covariance
            covSTrl = zeros(chans);
            for ti=1:length(spikelocs)
                tempdat = allTrlLFP{permtrl}(:,spikelocs(ti)-sets.nwin:spikelocs(ti)+sets.nwin); % take the EEG activity from all channels at timestamps surrounding the first spike
                
                spkLFPperm = [spkLFPperm; {tempdat}];
                
                tempdat = bsxfun(@minus,tempdat,mean(tempdat,2)); % demean data snippet
                covSTrl = covSTrl + (tempdat*tempdat')/sets.nwin; % calculate covariance of trough locked snippet and add it to existing covariance matrix
            end
            covS = covS + covSTrl; % average by number of added troughs
            
        end
        
        % this has to be done once a permutation
        covS  = covS  ./ numspks;
        
        [evecsP,evals] = eig(covS,bbcov);
        [~,compidx]    = sort(diag(evals)); % max component
        
        spkPowPerm = zeros(numspks, sniplen);
        for spk = 1 : size(spkLFPperm,1)
            spkLFPpermlin  = spkLFPperm{spk};
            spkLFPpermlin  = spkLFPpermlin' * evecsP(:,compidx(cmp));
            spkLFPpermlin  = [zeros(pddng,1); spkLFPpermlin; zeros(pddng,1)]; % padding
            spkPowPerm(spk,:)  = abs(fft(spkLFPpermlin)/sniplen);
        end
        
        trlPowPerm(perm,:) = mean(spkPowPerm,1);
%         disp(sprintf('Permutation #%d | elapsed time: %.1f seconds', perm, toc));
    end
    toc
    %% check if significant
    idx      = hz<=10;
    lowFperm = prctile(nanmean(trlPowPerm(:,idx),2),95, 1);
    lowF     = nanmean( nanmean(trlPow(:,idx),2) );
    
    idx      = and(hz>10, hz<=30);
    midFperm = prctile(nanmean(trlPowPerm(:,idx),2),95, 1);
    midF     = nanmean( nanmean(trlPow(:,idx),2) );
    
    idx      = and(hz>30, hz<=60);
    higFperm = prctile(nanmean(trlPowPerm(:,idx),2),95, 1);
    highF    = nanmean( nanmean(trlPow(:,idx),2) );
    
    if or(or(lowF >= lowFperm, midF >= midFperm), highF >= higFperm)
        signSFC = 1;
        evec    = evecs(:, compidx(cmp));
        comp = compidx(cmp);
    else
        signSFC = 0; % stays zero
        evec    = evecs(:, compidx(chans)); % take the biggest eigenvector
        comp = compidx(chans);
    end    
    
end

switch sets.server
    case 0
        savepath = ['X:\Luca\engram allocation lfp\SFC\method1\', num2str(filtf), filesep, num2str(sets.nwin), filesep];
    case 1
        savepath = ['/media/ldk898/rds-share/Luca/server/engram allocation lfp/', num2str(filtf), '/', num2str(sets.nwin), '/'];
end

save([savepath, 'SFC_', bidsID, '_', sesh, '_su', num2str(su)], 'trlPow', 'trlPowPerm', 'evecs', 'signSFC', 'evec', 'comp', '-v7.3')


%% Visualisation
% ______________

% vector of frequencies
% hz = linspace(0,srate,(sets.nwin*2)+1);

% figure
% for cmp = 1:chans
%
%     % plot signal
%     plotSig = nanmean(trlPow,1);
%     plot(hz, plotSig);
%     set(gca,'xlim',[0 90],'xtick',10:20:90), axis square
%     xlabel('Frequency (Hz)'), ylabel('Power (a.u.)')
%     title(sprintf('Comp #%d | #%d Spikes', cmp, numspks));
%
%     % plot 95% threshold
%     hold on
%     plotErr = prctile(trlPowPerm, 95,1);
%     plot(hz, plotErr);

% end

end % end of function