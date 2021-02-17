% params.split = 'no'; % 'yes'
% params.bl = 'wholeTrl'; % 'precue'
% params.meanDim = 'time'; % 'time'

function idxNdxAvgFiring(params)

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
dt = -1.375:0.001:10.375+0.001;
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
        spksSeg       = insertSpiketimes2(trig, clusterSpikes, [1 1], [-1.375 10.375])'; % 3 seconds prior to cue trigger until 1 second after response trigger
        
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
        
        
        %     switch sum(idxTrl)>0
        %         case 0 % ndx unit
        %             ndx = [ndx; mean(allConv,1)];
        %         case 1 % idx unit
        %             idxN = [idxN; mean(allConv(~idxTrl,:),1)];
        %             idxI = [idxI; mean(allConv( idxTrl,:),1)];
        %     end
        fprintf(repmat('\b',1,length(lineLength)-1));
        
    end
    
    switch period
        case 1 % encoding
            encConv = idxI;
        case 2 % retrieval
            retConv = idxI;
    end
    
   

    %% VISUALISATION
    % subplot(1,2,1); hold on;
    plotDT = [-1000:1:10000];
    % plotDat = nanmedian(ndx,1); %% SOME DO NOT FIRE DURING PRECUE PERIOD!!
    % plot(plotDT, plotDat, 'color', [0 0 0], 'linew', 2)
    %
    % shade     = std(ndx,0,1) / sqrt(size(ndx,1));
    shadeDT   = [plotDT, fliplr(plotDT)];
    % inBetween = [plotDat-shade, fliplr(plotDat+shade)];
    % fillHand = fill(shadeDT, inBetween, [0.1059, 0.6196, 0.4667] );
    % fillHand.FaceAlpha = 0.5;
    %
    %
    % xlim([-1000 10000])
    % ylim([-0.5 3.5])
    % xlabel('Time')
    % xticks([0 2000 5000 10000])
    % switch strcmp(encRet, 'enc')
    %     case 0
    %         xticklabels({'Cue', 'Q', '5s','10s'})
    %     case 1
    %         xticklabels({'Cue', 'Stim', '5s','10s'})
    % end
    % yticks(0)
    % ylabel('Firing rate (z-values)')
    % title('Non Index Units')
    % ax = gca;
    % ax.FontSize = 20;
    % ax.FontWeight = 'bold';
    % box off
    
    %%
    % IU INDEXED TRIALS
    subplot(1,2,period); hold on;
    plotDatI = mean(idxI,1);
    idxPlot = plot(plotDT, plotDatI, 'color', [0, 0, 0], 'linew', 2);
    
    % IU  NON-INDEXED TRIALS
    plotDatN = mean(idxN,1);
    ndxPlot = plot(plotDT, plotDatN, 'color', [0, 0, 0], 'linew', 2);
    
    % ERROR IDX
    shade     = std(idxI,0,1) / sqrt(size(idxI,1));
    inBetween = [plotDatI-shade, fliplr(plotDatI+shade)];
    fillHand = fill(shadeDT, inBetween, [0.4588, 0.4392, 0.7020] );
    fillHand.FaceAlpha = 0.5;
    
    % ERROR NDX
    shade     = std(idxN,0,1) / sqrt(size(idxN,1));
    inBetween = [plotDatN-shade, fliplr(plotDatN+shade)];
    fillHand = fill(shadeDT, inBetween, [0.1059, 0.6196, 0.4667] );
    fillHand.FaceAlpha = 0.5;
    
    % AXES
    xlim([-1000 10000])
    ylim([-0.5 5])
    xlabel('Time')
    xticks([0 2000 5000 10000])
    switch strcmp(encRet, 'enc')
        case 0
            xticklabels({'Cue', 'Q', '5s','10s'})
        case 1
            xticklabels({'Cue', 'Stim', '5s','10s'})
    end
    yticks([0:1:5])
    ylabel('Firing rate (z-values)')
    % title('Index Units')
    
%     % % LEGENDS
%     if strcmp(encRet, 'ret')
%         L(1) = plot(nan, nan, 'color', [0.4588, 0.4392, 0.7020]);
%         L(2) = plot(nan, nan, 'color', [0.1059, 0.6196, 0.4667]);
%         [legPos, hobj, ~, ~] = legend(L, {'Indexed Trials', 'Non-indexed Trials'}, 'FontSize',16, 'FontWeight', 'bold');
%         set(hobj,'LineWidth',15);
%         set(legPos, 'Position', [0.7629    0.7858    0.1962    0.1022])
%         legend('boxoff')
%     end
    
    % hobj = legend('a')
    % c=get(hobj,'Children')
    % set(c,'Color',[0 0 1])
    
    %     legend('boxoff');
    
    % switch strcmp(encRet, 'enc')
    %     case 0
    %         sgtitle('Retrieval', 'FontSize', 20, 'FontWeight', 'bold')
    %     case 1
    %         sgtitle('Encoding',  'FontSize', 20, 'FontWeight', 'bold')
    % end
    
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

saveTit = ['Split', params.split, '_BL', params.bl, '_MeanDim', params.meanDim, '.png'];
try
saveas(gcf, ['\\analyse4.psy.gla.ac.uk\project0309\Luca\visu\firingRate_idxNdx\', saveTit]);
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

%% XC

encWurst = reshape(encConv', 1, []);
retWurst = reshape(retConv', 1, []);
maxLag = 500;
[xc_emp, lags] = xcorr(encWurst, retWurst, maxLag, 'normalized');

% plot(lags, xc_emp, 'linew', 3, 'color', [0 0 0])
    
%% permtest
encRetConv = [encConv; retConv];
nperm = 500;
numTrl = size(encRetConv,1);
xc_perm = zeros(nperm, length(lags));

for perm = 1:nperm
    disp(perm);
    convPerm             = encRetConv( randperm( numTrl, numTrl), :);
    encPerm              = convPerm(1:numTrl/2,:);
    encPermWurst         = reshape(encPerm', 1, []);
    
    retPerm              = convPerm(numTrl/2+1:end, :);
    retPermWurst         = reshape(retPerm', 1, []);
    [xc_perm(perm,:), ~] = xcorr(encPermWurst, retPermWurst, maxLag, 'normalized');
end

meanPerm = mean(xc_perm,1);
stdPerm  = std(xc_perm, 0, 1);

xc_empNorm = xc_emp - meanPerm;
xc_empNorm = xc_empNorm ./ stdPerm;

figure(2); clf; hold on;
plot(lags, xc_empNorm, 'linew', 3, 'color', [0 0 0])
ax            = gca;
ax.FontSize   = 20;
ax.FontWeight = 'bold';
box off
title('Crosscorrelation between encoding and retrieval spiking (-1s : 10s)')
xlabel('Negative peak indicates earlier encoding peak')
    


end % END OF FUNCTION

% for ib = 1:2
%     if ib == 1
%         params.split = 'no';
%     else
%         params.split = 'yes';
%     end
%     
%     for ic = 1:2
%         if ic == 1
%             params.bl = 'wholeTrial';
%         else
%             params.bl = 'precue';
%         end
%         
%         
%         for id = 1:2
%             if id == 1
%                 params.meanDim = 'trl';
%             else
%                 params.meanDim = 'time';
%             end
%         end
%     end
%     
%     idxNdxAvgFiring(params)
%     
% end