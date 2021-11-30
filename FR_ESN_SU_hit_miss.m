function FR_ESN_SU_hit_miss

% gaussian kernel
mlength = 251;
mKernel = gausswin(mlength);

dt = -1.125:0.001:5.125+0.001;
dt = dt-0.0005; % center around 0

plotDT = [-1000:1:5000];
shadeDT   = [plotDT, fliplr(plotDT)];
suCount = 0;

%% fVSp
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\ESN_code_data\code_ESN\data\allSpks.mat', 'allSpks')
figure(1); %clf;
for cate = 1:3%:6
    for encRet = 1:2
        
        allConv = [];
        for su = 1:length(allSpks)
            suCount = suCount + 1;
            
            lineLength = fprintf('Cate %d/6 || encRet = %d || SU = %d out of %d (%d out of %d or %.2f%%).\n', cate, encRet, su, length(allSpks), suCount, length(allSpks)*2*6, (suCount/(length(allSpks)*2*6))*100 );
            
            switch cate
                case 1 % miss-ESN
                    
                    if allSpks(su).mESN == 0 % HAS TO BE A MISS ESN
                        continue
                    end
                    
                    himiDx                 = allSpks(su).himiDx; % 49 trials
                    reinstTrl              = allSpks(su).reinstTrl; % which of the miss trials are reinstated
                    idx                    = himiDx == 3 && reinstTrl;
                    
                    
                    trigEnc                = allSpks(su).encTrigger(idx,1);
                    trigRet                = allSpks(su).retTrigger(idx,1);
                    
                    
%                     if isempty(trigEnc)
%                         continue
%                     end
                    
                case 2 % HIT
                    
                    if allSpks(su).ESN == 0 %% HAS TO BE A HIT ESN
                        continue
                    end
                    
                    himiDx                 = allSpks(su).himiDx; % 49 trials
                    reinstTrl              = allSpks(su).reinstTrl; % which of the miss trials are reinstated
                    idx                    = himiDx == 1 && reinstTrl;
                    
                    trigEnc                = allSpks(su).encTrigger(idx,1);
                    trigRet                = allSpks(su).retTrigger(idx,1);

                case 3 % HALF HIT
                    
                    if allSpks(su).hmESN == 0 % HAS TO BE A MISS ESN
                        continue
                    end
                                        
                    himiDx                 = allSpks(su).himiDx; % 49 trials
                    reinstTrl              = allSpks(su).reinstTrl; % which of the miss trials are reinstated
                    idx                    = himiDx == 2 && reinstTrl;
                    
                    trigEnc                = allSpks(su).encTrigger(idx,1);
                    trigRet                = allSpks(su).retTrigger(idx,1);
                    
                    
                case 4 % hits SU
                    if allSpks(su).mESN == 1 || allSpks(su).ESN == 1 || allSpks(su).ESN == 2 || allSpks(su).hmESN == 1 % HAS TO BE A SU
                        continue
                    end
                    
                    if sum(allSpks(su).reinstTrl) > 0 || sum(allSpks(su).mReinstTrl) > 0 || sum(allSpks(su).hmReinstTrl) > 0
                        error('This should not happen');
                    end
                    
                    hitsIdx = allSpks(su).hitsIdx;
                    trigEnc = allSpks(su).encTrigger(hitsIdx,1);
                    trigRet = allSpks(su).retTrigger(hitsIdx,1);
                    
                case 5 % half miss SU
                    if allSpks(su).mESN == 1 || allSpks(su).ESN == 1 || allSpks(su).ESN == 2 || allSpks(su).hmESN == 1 % HAS TO BE A SU
                        continue
                    end
                    
                    missIdx = allSpks(su).himiDx;
                    missIdx = missIdx == 2;
                    trigEnc = allSpks(su).encTrigger(missIdx,1);
                    trigRet = allSpks(su).retTrigger(missIdx,1);
                    
                case 6 % (full) miss SU
                    
                    if allSpks(su).mESN == 1 || allSpks(su).ESN == 1 || allSpks(su).ESN == 2 || allSpks(su).hmESN == 1% HAS TO BE A SU
                        continue
                    end
                    
                    missIdx = allSpks(su).himiDx;
                    missIdx = missIdx == 3;
                    trigEnc = allSpks(su).encTrigger(missIdx,1);
                    trigRet = allSpks(su).retTrigger(missIdx,1);
                    
            end
            
            spkTms     = allSpks(su).spks/1000;
            
            if length(trigEnc) ~= length(trigRet)
                error('This should not happen');
            end
            
            if isempty(trigEnc) % this can happen if a patient has e.g. no misses
                continue
            end
            
            switch encRet
                case 1
                    trig = trigEnc;
                    mColor = [0.4588, 0.4392, 0.7020];
                case 2
                    trig = trigRet;
                    mColor = [0.1059, 0.6196, 0.4667];
            end
            
            spksSeg = insertSpiketimes2(trig, spkTms, [1 1], [-1.125 5.125])'; % 3 seconds prior to cue trigger until 1 second after response trigger
            
            if isempty(spksSeg)
                error('line 158')
            end
            
            subConv = [];
            for trl = 1 : length(trigEnc)
                
                [x,~]   = histcounts(spksSeg{trl}, dt);
                spkConv = conv(x, mKernel, 'same');
                
                subConv = [subConv; spkConv];
                
            end
            subConv = subConv(:,126:end-125); % cut off wings
            
            BL = subConv(:,1:1000); % precue period as baseline
            BL = mean(BL,2); % mean over time
            BLmean = mean(BL);
            BLstd  = std(BL, 0);
            BLstd = BLstd + 0.1;
            
            subConv = subConv - BLmean;
            subConv = subConv./BLstd;
            
            
            allConv = [allConv; subConv];
            
            fprintf(repmat('\b',1,length(lineLength)-1));
            
        end
        
        %% VISUALISATION
        subplot(3,2,cate); hold on;
        meanConv = mean(allConv,1);
        %         plot(plotDT, meanConv, 'linew', 3, 'color', mColor);
        plot(plotDT, meanConv, 'color', [0, 0, 0], 'linew', 2);
        
        
        % ERROR IDX
        shade     = std(allConv,0,1) / sqrt(size(allConv,1));
        inBetween = [meanConv-shade, fliplr(meanConv+shade)];
        fillHand  = fill(shadeDT, inBetween, mColor );
        fillHand.FaceAlpha = 0.5;
        
        
        
        %     % AXES
        %     ylim([-0.75 3.5])
        xlabel('Time')
        xticks([-1000 0 1000 2000 3000 5000 10000])
        xticklabels({'-1s', 'Cue', '1s', '2s', '3s', '5s','10s'})
        
        %     yticks([0:1:3])
        ylabel('Firing rate (z-values)')
        xlim([-1000 5000])
        
        trls = size(allConv,1);
        switch cate
            case 1
                title(sprintf('miss-ESNs (%d trials)', trls));
                
                if encRet == 1
                    output.missESN(1,:,:) = allConv;
                else
                    output.missESN(2,:,:) = allConv;
                end
                
            case 2
                title(sprintf('hit-ESNs (%d trials)', trls));
                if encRet == 1
                    output.hitESN(1,:,:) = allConv;
                else
                    output.hitESN(2,:,:)= allConv;
                end
                
            case 3
                title(sprintf('halfhit-ESNs (%d trials)', trls));
                if encRet == 1
                    output.himiESN(1,:,:) = allConv;
                else
                    output.himiESN(2,:,:) = allConv;
                end
                
            case 4
                title(sprintf('hit-SU (%d trials)', trls));
            case 5
                title(sprintf('halfhit-SU (%d trials)', trls));
            case 6
                title(sprintf('miss-SU (%d trials)', trls));
        end
        
    end % END OF encRet
end % END OF CATEGORIZATION

save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\outputESN.mat', 'output');

%% VISU
figure(1); clf;
for miss_himi = 1:2
    
    if miss_himi == 1 % hits vs. misses
        cc = 0;
    elseif miss_himi == 2 % hits vs himi
        cc = 2;
    end
    
    for encRet = 1:2
        
        subplot(2,2,encRet+cc) % ENCODING hitESN & missESN
        hold on;
        
        switch encRet
            case 1
                title('Encoding');
            case 2
                title('Retrieval')
        end
        
        dat1 = squeeze(output.hitESN (encRet,:,:)); % hitESNenc
        switch miss_himi
            case 1 % hits vs miss
                dat2 = squeeze(output.missESN(encRet,:,:)); % missESNenc
            case 2 % hits vs himi
                dat2 = squeeze(output.himiESN(encRet,:,:)); % missESNenc
        end
        
        for dats = 1:2
            switch dats
                case 1 % hits
                    dat = dat1;
                    mColor = [0.4588, 0.4392, 0.7020];
                case 2 % miss or himi
                    dat = dat2;
                    mColor = [0.1059, 0.6196, 0.4667];
            end
            
            meanConv = mean(dat,1);
            plot(plotDT, meanConv, 'color', [0, 0, 0], 'linew', 2);
            
            % ERROR IDX
            shade     = std(dat,0,1) / sqrt(size(dat,1));
            inBetween = [meanConv-shade, fliplr(meanConv+shade)];
            fillHand  = fill(shadeDT, inBetween, mColor );
            fillHand.FaceAlpha = 0.5;
            
    %% CLUSTER BASED STATISTICS
        clear cl1
    cl1.label = {'whatever'};
    for trl = 1:size(dat1,1)
        cl1.trial{trl} = dat1(trl,:);
        
        cl1.time(trl) = {-1000:1:5000};
        cl1.fample = 1000;
        
    end
    
    clear cl2
    cl2.label = {'whatever'};
    for trl = 1:size(dat2,1)
        cl2.trial{trl} = dat2(trl,:);
        
        cl2.time(trl) = {-1000:1:5000};
        cl2.fample = 1000;
        
    end
    
    cfg            = [];
    cfg.keeptrials = 'yes';
    timeLockCl1    = ft_timelockanalysis(cfg, cl1);
    timeLockCl2    = ft_timelockanalysis(cfg, cl2);
    
    cfg                  = [];
    cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
    cfg.statistic        = 'indepsamplesT'; % use the independent samples T-statistic as a measure to
    % evaluate the effect at the sample level
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;       % alpha level of the sample-specific test statistic that
    % will be used for thresholding
    cfg.clusterstatistic = 'maxsum';   % test statistic that will be evaluated under the
    % permutation distribution.
    cfg.tail             = 0;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025;      % alpha level of the permutation test
    cfg.numrandomization = 1000;       % number of draws from the permutation distribution
    cfg.design           = [ones(1,size(dat1,1)), ones(1,size(dat2,1))*2];
    [stat]               = ft_timelockstatistics(cfg, timeLockCl1, timeLockCl2);
    
    hstat                = stat.prob <= 0.05; % sign hstat as "1"
    tt                   = -1000:1:5000;
    tt(logical(hstat))   = 10;
    tt(~hstat)           = NaN; 
    
    plot([-1000:1:5000], tt, 'linew', 10, 'color', [1 0 0]);
    
    %% DEPENDENT SAMPLE
%     cstat = {};
%     
%     cstat{1}.label               = {'Channels'};
%     cstat{1}.time                = [-1.0:1/1000:5];
%     
%     cstat{1}.individual(:,1,:)   = dat1; %% WHY 3D?? IT'S ONLY (subj X time) classifier performance
%     cstat{1}.dimord              = 'subj_chan_time'; % 3D should be: 'subj_time_time'
%     % cstat{1}.avg                 = squeeze(mean(cstat{1}.individual)); % I THINK THIS IS THE AVERAGE OVER SUBJECTS
%     
%     cstat{2}                     = cstat{1};
%     cstat{2}.individual          = [];
%     cstat{2}.individual(:,1,:)   = dat2;
%     % cstat{2}.avg                 = squeeze(mean(cstat{2}.individual));
%     
%     
%     % FT stats
%     cfg                     = [];
%     cfg.latency             = [-1.0 5.0];
%     cfg.spmversion          = 'spm12';
%     cfg.channel             = 'all';
%     cfg.neighbours          = [];
%     cfg.minnbchan           = 0;
%     cfg.avgovertime         = 'no';
%     cfg.avgoverchan         = 'no';
%     cfg.computecritval      = 'yes';
%     cfg.parameter           = 'individual';
%     
%     cfg.statistic           = 'depsamplesT';
%     cfg.method              = 'montecarlo';
%     cfg.correctm            = 'cluster';
%     cfg.alpha               = .05;
%     cfg.clusteralpha        = .05;
%     cfg.correcttail         = 'no'; % one sided test
%     cfg.tail                = 1; % only test right tail
%     cfg.numrandomization    = 1000;
%     cfg.clusterstatistic    = 'maxsum';
%     cfg.clustertail         = cfg.tail;
%     
%     
%     % SO THIS SHOULD BE THE DESIGN MATRIX WHICH HAS TWO ENTRIES PER SUBJECT (EMP vs. PERM)
%     nSesh = size(dat1,1);
%     
%     % set up design matrix
%     cfg.design  = [1:nSesh, 1:nSesh; ones(1,nSesh), ones(1,nSesh)*2];
%     cfg.uvar    = 1;
%     cfg.ivar    = 2;
%     
%     % run stats
%     [stats] = ft_timelockstatistics(cfg, cstat{1}, cstat{2});
%     % stats.prob = 0 means difference
%     % stats.mask = 1 means difference
%     
%     % MARK SIGNIFICANT PERIODS
%     hstat                  = stats.mask == 1;
%     tvec                   = -1:1/sr:2;
%     
%     tvec(logical(hstat))   = 0.70;
%     tvec(~hstat)           = NaN;
%     
%     plot([-1:1/sr:2], tvec, 'linew', 10, 'color', [1 0 0]);
%     
            
        end
    end
    
    % LEGEND
    L(1) = plot(nan, nan, 'color', [0.4588 0.4392 0.7020]);
    L(2) = plot(nan, nan, 'color', [0.1059 0.6196 0.4667]);
    switch miss_himi
        case 1 % miss
            
            [legPos, hobj, ~, ~] = legend(L, {'Hits', 'Miss'}, 'FontSize',15, 'FontWeight', 'bold');
            set(hobj,'LineWidth',15);
            %             set(legPos, 'Position', [0.7629    0.7858    0.1962    0.1022])
            legend('boxoff')
            
        case 2
            
            [legPos, hobj, ~, ~] = legend(L, {'Hits', 'Half Hits'}, 'FontSize',15, 'FontWeight', 'bold');
            set(hobj,'LineWidth',15);
            legend('boxoff')
              
    end
end



%% LIKE FOR ESNs
clear
folderpath = '\\analyse4.psy.gla.ac.uk\project0309\Luca\ESN_code_data\code_ESN\data'; % wherever you have downloaded this folder
addpath(genpath(folderpath))
load([folderpath, '\allSpks.mat'], 'allSpks');
outputH = struct; outputM = struct;
counterH = 0; counterM = 0;
tic
for himi = 1:2
    for spk = 1 : size(allSpks,2)
        timePassed = toc;
        lineLength = fprintf('Analysing SU#%d of %d | Time passed: %.1f min \n', spk, size(allSpks,2), (timePassed/60));
        % GET THE SPIKES
        spkTms = allSpks(spk).spks;
        
        % ENC: CUE-RESP || RET: CUE-RESP
        encTrigger = round(allSpks(spk).encTrigger(:,[1 3])*1000);     % CUE ONSET UNTIL RESPONSE
        retTrigger = round(allSpks(spk).retTrigger(:,[1])*1000);       % CUE ONSET UNTIL...
        retTrigger = [ retTrigger round(allSpks(spk).retRT(:)*1000) ]; % ... RESPONSE
        
        % EXCLUDE CN PERIOD
        encTrigger(:,1) = encTrigger(:,1) + 2000; % NOW STARTING AT RESPONSE ONSET
        switch himi
            case 1 % hits
                himiDx     = allSpks(spk).himiDx == 1; % hits
                reinstTrl  = allSpks(spk).reinstTrl';
                idx        = himiDx & reinstTrl;
                
            case 2 % full misses
                himiDx   = allSpks(spk).himiDx == 3; % full misses
                reinstTrl  = allSpks(spk).reinstTrl';
                idx        = himiDx & reinstTrl;
        end
        
        % SEGMENT THE SPIKETIMES INTO TRIALS
        spkTrlEnc    = [];
        spkTrlRet    = [];
        for trl = 1:size(encTrigger,1)
            % ENCODING
            trlLen         = (encTrigger(trl,2) - encTrigger(trl,1)) / 1000;
            spkTrl         = spkTms(spkTms>=encTrigger(trl,1) & spkTms<encTrigger(trl,2)) - encTrigger(trl,1) + 1;
            spkTrl         = size(spkTrl,1);
            spkTrlEnc(trl) = spkTrl / trlLen; % TRANSFORM IN Hz
            % RETRIEVAL
            trlLen         = (retTrigger(trl,2) - retTrigger(trl,1))  / 1000;
            spkTrl         = spkTms(spkTms>=retTrigger(trl,1) & spkTms<retTrigger(trl,2)) - retTrigger(trl,1) + 1;
            spkTrl         = size(spkTrl,1);
            spkTrlRet(trl) = spkTrl / trlLen;
        end
        
        % Z-SCORE SPIKE NUMBER
        spkTrlEnc = (spkTrlEnc - mean(spkTrlEnc)) / std(spkTrlEnc);
        spkTrlRet = (spkTrlRet - mean(spkTrlRet)) / std(spkTrlRet);
        ewp       = spkTrlEnc .* spkTrlRet; % ELEMENT WISE PRODUCT
        
        switch himi
            case 1
                if allSpks(spk).ESN == 1
                    counterH = counterH + 1;
                    outputH(counterH).ESNenc = spkTrlEnc(idx);
                    outputH(counterH).ESNret = spkTrlRet(idx);
                    outputH(counterH).ESNewp = ewp(idx);
                end
            case 2
                if allSpks(spk).mESN == 1
                    counterM = counterM + 1;
                    outputM(counterM).ESNenc = spkTrlEnc(idx);
                    outputM(counterM).ESNret = spkTrlRet(idx);
                    outputM(counterM).ESNewp = ewp(idx);
                end
        end
    end
end

figure(1); clf;
bar([outputH.ESNenc])
bar(mean([outputH.ESNenc]))

figure(1); clf;
bar([mean([outputH.ESNenc]) mean([outputH.ESNret]) mean([outputM.ESNenc]) mean([outputM.ESNret]) ])
xticklabels({'ESN Encoding', 'ESN Retrieval', 'miss-ESN Encoding', 'missESN Retrieval'})
figure(2); clf;
bar([mean([outputH.ESNewp]) mean([outputM.ESNewp])])
xticklabels({'ESN EWP', 'miss-ESN EWP'})
ttest([outputH.ESNenc], [outputM.ESNenc])
ttest2([outputH.ESNenc], [outputM.ESNenc])

[h,p,~]=ttest2([outputH.ESNenc], [outputM.ESNenc])
[h,p,~]=ttest2([outputH.ESNewp], [outputM.ESNewp])
[h,p,~]=ttest2([outputH.ESNret], [outputM.ESNret])


end % END OF FUNCTION