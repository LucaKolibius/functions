function idxNdxAvgFiring_vs_CN
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\declan tuning\allTun.mat', 'allTun'); % created with extractTunings
plotTun = mean(allTun,1);
shadeTun   = std(allTun,0,1) / sqrt(size(allTun,1));

params.split = 'no'; % 'yes'
params.bl = 'wholeTrl'; % 'precue'
params.meanDim = 'time'; % 'time'

%% SETUP FIGURE
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
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')

% encRet = params.encRet;
% dt = -1:0.001:10;
dt = -0.375:0.001:1.375+0.001;
dt = dt-0.0005; % center around 0


%% gaussian kernel
mlength = [-0.075:0.0002:0.075];
mSigma  = 0.02; % 20ms
% mSigma  = 0.01;
mKernel = normpdf(mlength,0,mSigma);
mKernel = mKernel/max(mKernel); % normalize peak to 1

for period = 1:2
    if period == 1
        encRet = 'enc';
    else
        encRet = 'ret';
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
        spksSeg       = insertSpiketimes2(trig, clusterSpikes, [1 1], [-0.375 1.375])'; % 3 seconds prior to cue trigger until 1 second after response trigger
        
        allConv = [];
        for trl = 1 : length(trig)
            
            [x,~]                = histcounts(spksSeg{trl}, dt);
            spkConv              = conv(mKernel, x);
            
            spkConv(1:750)       = [];
            spkConv(end-749:end) = [];
            
            allConv = [allConv; spkConv];
            
        end
        
        switch strcmp(params.bl, 'precue')
            case 0 %% use whole trial for baseline
                BL = allConv;
            case 1 %% use precue period for baseline
                BL = allConv(:,1:1000);
        end
        
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
                allConvIdx = allConv( idxTrl,:);
                allConvNdx = allConv(~idxTrl,:);
                
                BLidx = BL( idxTrl,:);
                BLndx = BL(~idxTrl,:);
                
                switch strcmp(params.meanDim, 'trl')
                    case 0 % over time
                        BLidx = mean(BLidx,2);
                        BLndx = mean(BLndx,2);
                    case 1 % over trl
                        BLidx = mean(BLidx,1);
                        BLndx = mean(BLndx,1);
                end
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
                BLmeanIdx = mean(BLidx);
                BLstdIdx  = std( BLidx, 0);
                
                BLmeanNdx = mean(BLndx);
                BLstdNdx  = std( BLndx, 0);
                
                if strcmp(params.bl, 'precue')
                    BLstdIdx = BLstdIdx + 0.1;
                    BLstdNdx = BLstdNdx + 0.1;
                end
                
                allConvIdx = allConvIdx - BLmeanIdx;
                allConvIdx = allConvIdx ./ BLstdIdx;
                
                allConvNdx = allConvNdx - BLmeanNdx;
                allConvNdx = allConvNdx ./ BLstdNdx;
        end
        
        idxN = [idxN; allConvNdx];
        idxI = [idxI; allConvIdx];
        
        
    end

    %% VISUALISATION
    plotDT  = [0:1:1000];
    shadeDT = [plotDT, fliplr(plotDT)];
    subplot(1,2,period); hold on;
    
    % CONCEPT UNITS  
    plot(plotDT, plotTun, '--', 'color', [0 0 0], 'linew', 2) 
    inBetween = [plotTun-shadeTun, fliplr(plotTun+shadeTun)];
    fillHand  = fill(shadeDT, inBetween, [0 0 0]);
    fillHand.FaceAlpha = 0.5;
    
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
    fillHand.FaceAlpha = 0.5;
    
    % ERROR NDX
    shade     = nanstd(idxN,0,1) / sqrt(size(idxN,1));
    inBetween = [plotDatN-shade, fliplr(plotDatN+shade)];
    fillHand  = fill(shadeDT, inBetween, [0.1059, 0.6196, 0.4667] );
    fillHand.FaceAlpha = 0.5;
    
    % AXES
    xlim([-0 1000])
    ylim([-0.5 2])
    xlabel('Time [ms]')
    xticks([0 500 1000])
    xticklabels({'Cue' '500' '1000' })
%     yticks([0:0.1:1])
    ylabel('Firing rate (z-values)')
    
    %% T-TEST
    pTT = [];
    for tp = 1:1001
        [pTT(tp),~] = ttest2(allTun(:,tp), idxI(:,tp));
    end
    tt = 1:1:1001;
    tt(~pTT) = NaN;
    tt(logical(pTT))  = 1.5;
    
    plot([1:1:1001], tt, 'linew', 3, 'color', [1 0 0]);

    %%
    mtit = ['Split into idx/ndx: ', params.split, ' | Baseline: ', params.bl, ' | MeanDim: ', params.meanDim];
    sgtitle(mtit);

    switch strcmp(encRet, 'enc')
        case 0
            title('Retrieval')
        case 1
            title('Encoding')
    end
    
    ax            = gca;
    ax.FontSize   = 20;
    ax.FontWeight = 'bold';
    box off
    
end % END FOR ENC / RET PERIOD

L(1) = plot(nan, nan, 'color', [0.4588, 0.4392, 0.7020]);
L(2) = plot(nan, nan, 'color', [0.1059, 0.6196, 0.4667]);
L(3) = plot(nan, nan, 'color', [0 0 0]);
[legPos, hobj, ~, ~] = legend(L, {'Indexed Trials', 'Non-indexed Trials', 'Concept units'}, 'FontSize',16, 'FontWeight', 'bold');
set(hobj,'LineWidth',15);
set(legPos, 'Position', [0.7629    0.7858    0.1962    0.1022])
legend('boxoff')

saveTit = ['Split', params.split, '_BL', params.bl, '_MeanDim', params.meanDim, '.png'];
try
saveas(gcf, ['\\analyse4.psy.gla.ac.uk\project0309\Luca\visu\firingRate_idxNdx_vs_CN\', saveTit]);
catch
    disp('cannot save image')
end

% %% idxI seems to be delta modulated!!
%
% tt= mean(idxI,1);
% pow = abs(fft(tt));
%
% fftDT = linspace(0,1000, length(tt));
% plot(fftDT, pow);
% xlim([1 20])
%
% plot(dt(2:end), tt)




end % END OF FUNCTION
