function idxNdxAvgFiring_vs_CN
addpath(' \\analyse4.psy.gla.ac.uk\project0309\Luca') % sendmail is in here
colEnc         = [0.368627450980392, 0.235294117647059, 0.60000000000000000];
colRet         = [0.901960784313726, 0.380392156862745, 0.00392156862745098];
colEncNdx      = [0.929411764705882, 0.972549019607843, 0.98431372549019600];
colRetNdx      = [0.996078431372549, 0.941176470588235, 0.85098039215686300];

%% CORRECT THE RIGHT TUNING DATA!!
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\declan tuning\allTun_blPreCue.mat', 'allTun'); % created with extractTunings
plotTun        = mean(allTun,1);
shadeTun       = std(allTun,0,1) / sqrt(size(allTun,1));

params.split   = 'no'; % 'yes'
params.bl      = 'precue'; % 'wholeTrl'
params.meanDim = 'time'; % 'time'

% %% SETUP FIGURE
mFigH = figure(1); clf;
%     set(gcf, 'units','normalized','outerposition',[0 0 1 1])

MP = get(0, 'MonitorPositions');
N = size(MP, 1);
% Might want to set an initial position this to some reasonable location
% in the event of the window being "Restored Down".
newPosition = MP(1,:);

if size(MP, 1) == 1
    % Single monitor
    set(mFigH,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
else
    % Multiple monitors - shift to the Nth monitor.
    newPosition(1) = newPosition(1) + MP(N,1);
end
mFigH.set('Position', newPosition, 'units', 'normalized');
mFigH.WindowState = 'maximized'; % Maximize with respect to current monitor.

%% LOAD IN DATA AND SET PARAMETER
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp.mat', 'allSpks')

% encRet = params.encRet;
% dt = -1:0.001:10;
dt = -0.625:0.001:1.125+0.001;
dt = dt-0.0005; % center around 0


%% gaussian kernel
% mlength = [-0.075:0.0002:0.075];
% mSigma  = 0.02; % 20ms
% % mSigma  = 0.01;
% mKernel = normpdf(mlength,0,mSigma);
% mKernel = mKernel/max(mKernel); % normalize peak to 1
mlength  = 251;
mKernel  = gausswin(mlength);
saveProb = [];

for period = 1:2
%     %% SETUP FIGURE
%     mFigH = figure(period); clf;
%     %     set(gcf, 'units','normalized','outerposition',[0 0 1 1])
%     
%     MP = get(0, 'MonitorPositions');
%     N = size(MP, 1);
%     % Might want to set an initial position this to some reasonable location
%     % in the event of the window being "Restored Down".
%     newPosition = MP(1,:);
%     
%     if size(MP, 1) == 1
%         % Single monitor
%         set(mFigH,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
%     else
%         % Multiple monitors - shift to the Nth monitor.
%         newPosition(1) = newPosition(1) + MP(N,1);
%     end
%     mFigH.set('Position', newPosition, 'units', 'normalized');
%     mFigH.WindowState = 'maximized'; % Maximize with respect to current monitor.
    
    if period == 1
        encRet    = 'enc';
        curCol    = colEnc;
        curColNdx = colEncNdx;
    else
        encRet    = 'ret';
        curCol    = colRet;
        curColNdx = colRetNdx;
    end
    
    
    ndx  = [];
    idxN = [];
    idxI = [];
    for su = 1: length(allSpks)
        %% ignore SU for now
        if allSpks(su).iu == 0
            continue
        end
        
        lineLength = fprintf('%d SU out of %d SU done (%.2f%%).\n', su, length(allSpks), su/length(allSpks)*100);
        
        %% TRIGGER IS IN SECONDS
        switch strcmp(encRet, 'enc')
            case 0 % RETRIEVAL TRIGGER
                trig = allSpks(su).retTrigger(allSpks(su).hitsIdx, [1]);
            case 1 % ENCODING TRIGGER
                trig = allSpks(su).encTrigger(allSpks(su).hitsIdx, [1]);
        end
        
        idxTrl     = allSpks(su).idxTrl;
        
        %% SPIKES IN SECONDS
        clusterSpikes = allSpks(su).spks/1000;
        spksSeg       = insertSpiketimes2(trig, clusterSpikes, [1 1], [-0.625 1.125])'; % 3 seconds prior to cue trigger until 1 second after response trigger
        
        allConv = [];
        for trl = 1 : length(trig)
            
            [x,~]                = histcounts(spksSeg{trl}, dt);
            spkConv              = conv(x, mKernel, 'same');
            
            %             spkConv(1:750)       = [];
            %             spkConv(end-749:end) = [];
            
            allConv = [allConv; spkConv];
            
        end
        allConv = allConv(:,126:end-125); % cut off bleeding from taking bigger TOI (to prevent edge effects)
        
        switch strcmp(params.bl, 'precue')
            case 0 %% use whole trial for baseline
                BL = allConv(:,501:end);
            case 1 %% use precue period for baseline
                BL = allConv(:,1:500);
        end
        allConv = allConv(:,501:end);

        switch strcmp(params.split, 'yes')
            case 0 % don't split into idx and ndx
                %             BL = BL;
                
                switch strcmp(params.meanDim, 'trl')
                    case 0
                        BL = mean(BL,2);
                    case 1
                        BL = mean(BL,1);
                end
                
            case 1 % split into idx/ndx
%                 allConvIdx = allConv( idxTrl,:);
%                 allConvNdx = allConv(~idxTrl,:);
%                 
%                 BLidx = BL( idxTrl,:);
%                 BLndx = BL(~idxTrl,:);
%                 
%                 switch strcmp(params.meanDim, 'trl')
%                     case 0 % over time
%                         BLidx = mean(BLidx,2);
%                         BLndx = mean(BLndx,2);
%                     case 1 % over trl
%                         BLidx = mean(BLidx,1);
%                         BLndx = mean(BLndx,1);
%                 end
        end
        
        
        switch strcmp(params.split, 'yes')
            case 0 % no split
                BLmean = mean(BL);
                BLstd  = std(BL, 0);
                
                if strcmp(params.bl, 'precue')
                    BLstd = BLstd + 0.1;
                end
                
                allConv = allConv - BLmean;
                allConv = allConv./BLstd;
                
                allConvIdx = allConv( idxTrl,:);
                allConvNdx = allConv(~idxTrl,:);
            case 1 % split into idx/ndx
%                 BLmeanIdx = mean(BLidx);
%                 BLstdIdx  = std( BLidx, 0);
%                 
%                 BLmeanNdx = mean(BLndx);
%                 BLstdNdx  = std( BLndx, 0);
%                 
%                 if strcmp(params.bl, 'precue')
%                     BLstdIdx = BLstdIdx + 0.1;
%                     BLstdNdx = BLstdNdx + 0.1;
%                 end
%                 
%                 allConvIdx = allConvIdx - BLmeanIdx;
%                 allConvIdx = allConvIdx ./ BLstdIdx;
%                 
%                 allConvNdx = allConvNdx - BLmeanNdx;
%                 allConvNdx = allConvNdx ./ BLstdNdx;
        end
        
        idxN = [idxN; allConvNdx];
        idxI = [idxI; allConvIdx];
        
        fprintf(repmat('\b', 1, lineLength));
    end
    
    %% VISUALISATION
    plotDT  = [0:1:1000];
    shadeDT = [plotDT, fliplr(plotDT)];
        subplot(1,2,period);
    hold on;
    
    % CONCEPT UNITS
    plot(plotDT, plotTun, '--', 'color', [0 0 0], 'linew', 2)
    inBetween = [plotTun-shadeTun, fliplr(plotTun+shadeTun)];
    fillHand  = fill(shadeDT, inBetween, [0 0 0]);
    fillHand.FaceAlpha = 0.25;
    
    % IU INDEXED TRIALS
    plotDatI = nanmean(idxI,1);
    idxPlot  = plot(plotDT, plotDatI, 'color', [0, 0, 0], 'linew', 2);
    
    % IU  NON-INDEXED TRIALS
    plotDatN = nanmean(idxN,1);
    ndxPlot  = plot(plotDT, plotDatN, 'color', [0, 0, 0], 'linew', 2);
    
    % ERROR IDX
    shade     = nanstd(idxI,0,1) / sqrt(size(idxI,1));
    inBetween = [plotDatI-shade, fliplr(plotDatI+shade)];
    fillHand  = fill(shadeDT, inBetween, [0.4588, 0.4392, 0.7020] );
    %     fillHand  = fill(shadeDT, inBetween, curCol ); % curCol = colEnc or colRet
    fillHand.FaceAlpha = 0.5;
    
    % ERROR NDX
    shade     = nanstd(idxN,0,1) / sqrt(size(idxN,1));
    inBetween = [plotDatN-shade, fliplr(plotDatN+shade)];
    fillHand  = fill(shadeDT, inBetween, [0.1059, 0.6196, 0.4667] );
    %     fillHand  = fill(shadeDT, inBetween, curColNdx);
    fillHand.FaceAlpha = 0.5;
    
    % AXES
    xlim([-0 1000])
    ylim([-0.5 2])
    xlabel('Time [ms]')
    xticks([0 500 1000])
    xticklabels({'Cue' '500' '1000' })
    %     yticks([0:0.1:1])
    ylabel('Firing rate (z-values)')
    
    %% RMANOVA (would need binning)
    %     grp = [ones(size(allTun,1), 1); zeros(size(idxI,1),1)];
    %     mTab = table(grp, 'VariableNames', {'grp'});
    %
    %     allDat = [allTun; idxI];
    %     for tp = 1:size(allTun,2)
    %         newTab = table(allDat(:,tp), 'VariableNames', {sprintf('tp%d', tp)});
    %         mTab   = [mTab newTab];
    %     end
    %     tdt = 1:1:1001;
    %     rm = fitrm(mTab,'tp1-tp1001 ~ grp','WithinDesign',tdt);
    %     ranovatbl = ranova(rm)
    %
    %     for ii = 1:5
    %     disp(period)
    %     disp(ranovatbl)
    %     end
    
    %% STATISTICS: cluster based test
% %     clear tun
% %     tun.label = {'whatever'};
% %     for trl = 1:size(allTun,1)
% %         tun.trial{trl} = allTun(trl,:);
% %         
% %         tun.time(trl) = {0:1:1000};
% %         tun.fample = 1000;
% %         
% %     end
% %     
% %     clear idx
% %     idx.label = {'whatever'};
% %     for trl = 1:size(idxI,1)
% %         idx.trial{trl} = idxI(trl,:);
% %         
% %         idx.time(trl) = {0:1:1000};
% %         idx.fample = 1000;
% %         
% %     end
% %     
% %     cfg            = [];
% %     cfg.keeptrials = 'yes';
% %     timeLockTun    = ft_timelockanalysis(cfg, tun);
% %     timeLockIdx    = ft_timelockanalysis(cfg, idx);
% %     
% %     cfg                  = [];
% %     cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
% %     cfg.statistic        = 'indepsamplesT'; % use the independent samples T-statistic as a measure to
% %     % evaluate the effect at the sample level
% %     cfg.correctm         = 'cluster';
% %     cfg.clusteralpha     = 0.05;       % alpha level of the sample-specific test statistic that
% %     % will be used for thresholding
% %     cfg.clusterstatistic = 'maxsum';   % test statistic that will be evaluated under the
% %     % permutation distribution.
% %     cfg.tail             = 0;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
% %     cfg.clustertail      = 0;
% %     cfg.alpha            = 0.025;      % alpha level of the permutation test
% %     cfg.numrandomization = 10000;       % number of draws from the permutation distribution
% %     cfg.design           = [ones(1,size(allTun,1)), ones(1,size(idxI,1))*2];
% %     [stat]               = ft_timelockstatistics(cfg, timeLockTun, timeLockIdx);
% %     
% %     hstat                = stat.prob <= 0.05; % sign hstat as "1"
% %     tt                   = 0:1:1000;
% %     
% %     saveProb = [saveProb; stat.prob];
% %     switch period
% %         case 1 % enc
% %             tt(logical(hstat))   = 2.5;
% %             tt(~hstat) = NaN;
% %         case 2 % ret
% %             tt(logical(hstat))   = 2.5;
% %             tt(~hstat) = NaN;
% %             
% %     end
% %     
% %     plot([0:1:1000], tt, 'linew', 4, 'color', [1 0 0]);
% %     ylim([-0.5 3])
    %     figure(10); hold on; subplot(2,1,period); imagesc(stat.prob)
    %% T-TEST
    %     hTT = [];
    %     for tp = 1:1001
    %         [hTT(tp),~,~,stats] = ttest2(allTun(:,tp), idxI(:,tp));
    %         tTT(tp) = abs(stats.tstat);
    %     end
    %     tt       = 1:1:1001;
    %     tt(~hTT) = NaN;
    %     tt(logical(hTT))  = 1.5;
    %
    %     plot([1:1:1001], tt, 'linew', 3, 'color', [1 0 0]);
    
    %%
    mtit = ['Split into idx/ndx: ', params.split, ' | Baseline: ', params.bl, ' | MeanDim: ', params.meanDim];
    %     sgtitle(mtit);
    
    switch strcmp(encRet, 'enc')
        case 0
            title('Retrieval')
            L(1) = plot(nan, nan, 'color', [0.4588 0.4392 0.7020]);
            L(2) = plot(nan, nan, 'color', [0.1059 0.6196 0.4667]);
            L(3) = plot(nan, nan, 'color', [0 0 0]);
            [legPos, hobj, ~, ~] = legend(L, {'Reinstated Trials', 'Non-reinstated Trials', 'Tuned units'}, 'FontSize',25, 'FontWeight', 'bold');
            set(hobj,'LineWidth',15);
            set(legPos, 'Position', [0.7629    0.7858    0.1962    0.1022])
            legend('boxoff')
        case 1
            title('Encoding')
    end
    
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',25,'FontWeight','Bold', 'LineWidth', 2);
    box off
    
    
   
    
%     save(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\idxNdxFireCN_',encRet,'.mat'],'-v7.3');

    
end % END FOR ENC / RET PERIOD

    saveTit = ['Split', params.split, '_BL', params.bl, '_MeanDim', params.meanDim, '.svg'];
    try
        saveas(gcf, ['\\analyse4.psy.gla.ac.uk\project0309\Luca\visu\firingRate_idxNdx_vs_CN\', saveTit]);
    catch
        disp('cannot save image')
    end

    
% save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\pVal_firingRate_idxVScn', 'saveProb')
sendEmail('Script firing rate done', 'Firing Rate IDX/CN done for ENC & RET. 10000perm.')

end % END OF FUNCTION
