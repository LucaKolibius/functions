function idxNdxAvgFiring_XC

params.split   = 'no'; % 'yes'
params.bl      = 'wholeTrl'; % 'precue'
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
% dtold = 0:0.001:3;
dt = -0.375:0.001:3.376;% !!!! Changed !!!!!!!
dt = dt-0.0005; % center around 0


%% gaussian kernel
%mlength = [-0.125:0.001:0.125];% !!!!! changed !!!!! the Kernel seems to much too long, i.e. it appears to be 750 ms!! I would normally not expect this to be longer than +/- 125 ms
%mSigma  = 0.02; % 20ms
% mSigma  = 0.01;
%mKernel = normpdf(mlength,0,mSigma);
%mKernel = mKernel/max(mKernel); % normalize peak to 1
mlength = 251;
mKernel = gausswin(mlength);
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
        spksSeg       = insertSpiketimes2(trig, clusterSpikes, [1 1], [-0.375 3.375])'; % -0.375 seconds prior to cue trigger until 3.375 seconds after cue trigger
        
        
        %% CONVOLVE SPIKES
        allConv = [];
        for trl = 1 : length(trig)
            
            [x,~]                = histcounts(spksSeg{trl}, dt);
            %spkConv              = conv(mKernel, x);
            spkConv              = conv(x,mKernel,'same'); % same prevents kernel bleeding
%            spkConv(1:floor(numel(mKernel)/2)-1)       = []; 
%            spkConv(end-floor(numel(mKernel)/2):end) = [];
            
            allConv = [allConv; spkConv];
            
        end
        allConv = allConv(:,376:end-375); % cut off bleeding from taking bigger TOI (to prevent edge effects)
        
        %% NORMALISATION %% I don't understand why you normalize? I don't think you need to do this! I therefore switched it off
        switch strcmp(params.bl, 'precue')
            case 0 %% use whole trial for baseline
                BL = allConv;
            case 1 %% use precue period for baseline
                BL = allConv(:,1:1000);
        end
        
        switch strcmp(params.split, 'yes')
            case 0 % don't split into idx and ndx  % <-
                %             BL = BL;
                
                switch strcmp(params.meanDim, 'trl')
                    case 0
                        BL = mean(BL,2); % <-
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
            case 0 % no split % <-
                BLmean = mean(BL);
                BLstd  = std(BL, 0);
                
                if strcmp(params.bl, 'precue')
                    BLstd = BLstd + 0.1;
                end
                
                % simons approach:
                % allConv = allConv;% !!!!! changed !!!!!
                % allConv = allConv;% !!!!! changed !!!!! ./BLstd;
                allConv = allConv -  BLmean;
                allConv = allConv ./ BLstd;
                
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
                
                allConvIdx = allConvIdx;% !!!!! changed !!!!! - BLmeanIdx;
                allConvIdx = allConvIdx;% !!!!! changed !!!!! ./ BLstdIdx;
                
                allConvNdx = allConvNdx;% !!!!! changed !!!!! - BLmeanNdx;
                allConvNdx = allConvNdx;% !!!!! changed !!!!! ./ BLstdNdx;
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
        
    end % END OF SU LOOP
    
    %% STORE CONVOLUTION
    switch period
        case 1 % encoding
            encConv = idxN;
        case 2 % retrieval
            retConv = idxN;
    end
    
    %% VISUALISATION
    plotDT    = [0:1:3000];
    shadeDT   = [plotDT, fliplr(plotDT)];

    % IU INDEXED TRIALS
    subplot(1,2,period); hold on;
    plotDatI = mean(idxI,1);
    idxPlot  = plot(plotDT, plotDatI, 'color', [0, 0, 0], 'linew', 2);
    
    % IU  NON-INDEXED TRIALS
    plotDatN = mean(idxN,1);
    ndxPlot  = plot(plotDT, plotDatN, 'color', [0, 0, 0], 'linew', 2);
    
    % ERROR IDX
    shade     = std(idxI,0,1) / sqrt(size(idxI,1));
    inBetween = [plotDatI-shade, fliplr(plotDatI+shade)];
    fillHand  = fill(shadeDT, inBetween, [0.4588, 0.4392, 0.7020] );
    fillHand.FaceAlpha = 0.5;
    
    % ERROR NDX
    shade     = std(idxN,0,1) / sqrt(size(idxN,1));
    inBetween = [plotDatN-shade, fliplr(plotDatN+shade)];
    fillHand  = fill(shadeDT, inBetween, [0.1059, 0.6196, 0.4667] );
    fillHand.FaceAlpha = 0.5;
    
    % AXES
    xlim([-0 3000])
    ylim([-0.5 3])
    xlabel('Time')
    xticks([0 2000 3000])
    switch strcmp(encRet, 'enc')
        case 0
            xticklabels({'Cue', 'Q', '3s'})
        case 1
            xticklabels({'Cue', 'Stim', '3s'})
    end
    yticks([0:1:5])
    ylabel('Firing rate (z-values)')
    
%     % % LEGENDS
%     if strcmp(encRet, 'ret')
%         L(1) = plot(nan, nan, 'color', [0.4588, 0.4392, 0.7020]);
%         L(2) = plot(nan, nan, 'color', [0.1059, 0.6196, 0.4667]);
%         [legPos, hobj, ~, ~] = legend(L, {'Indexed Trials', 'Non-indexed Trials'}, 'FontSize',16, 'FontWeight', 'bold');
%         set(hobj,'LineWidth',15);
%         set(legPos, 'Position', [0.7629    0.7858    0.1962    0.1022])
%         legend('boxoff')
%     end
    
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
% 
% %% SAVE IMAGE
% saveTit = ['Split', params.split, '_BL', params.bl, '_MeanDim', params.meanDim, '.png'];
% try
% saveas(gcf, ['\\analyse4.psy.gla.ac.uk\project0309\Luca\visu\firingRate_idxNdx_shortTrial\', saveTit]);
% catch
%     disp('cannot save image')
% end


% %% XC
% 
% encWurst = reshape(encConv', 1, []);
% retWurst = reshape(retConv', 1, []);
% maxLag   = 500;
% [xc_emp, lags] = xcorr(encWurst, retWurst, maxLag);
% [~,eloc] = findpeaks(xc_emp, 'SortStr','descend', 'NPeaks', 1);
% eloc=lags(eloc);
% % plot(lags, xc_emp, 'linew', 3, 'color', [0 0 0])
%     
% %% PERMTEST
% encPerm=encConv;
% retPerm=retConv;
% nperm      = 1000;
% numTrl     = size(encPerm,1);
% xc_perm    = zeros(nperm, length(lags));
% 
% for perm = 1:nperm
%     flp=sign(rand(1,numTrl)-0.5);
%     disp(perm);
%     for k=1:numTrl
%         % randomly swap encoding with retrieval but preserve trial
%         % structure
%         if flp(k)==1
%             encPerm(k,:) = encConv(k,:);
%             retPerm(k,:) = retConv(k,:);
%         else
%             encPerm(k,:) = retConv(k,:);
%             retPerm(k,:) = encConv(k,:);
%         end
%     end
%     encPermWurst         = reshape(encPerm', 1, []);
%     retPermWurst         = reshape(retPerm', 1, []);
%     [xc_perm(perm,:), ~] = xcorr(encPermWurst, retPermWurst, maxLag);
%     [rpk,rloc] = findpeaks(xc_perm(perm,:), 'SortStr','descend', 'NPeaks', 1);
%     rlocs(perm,1)=lags(rloc);
% end
% 
% %% NORMALIZE XC WITH PERMUTATION XC
% meanPerm   = mean(xc_perm,1);
% meanPerm2   =mean(xc_perm,2); % !!!! Average across columns to get std across lags across randomization runs 
% stdPerm    = std(meanPerm2, 0, 1);
% 
% xc_empNorm = xc_emp - meanPerm;
% xc_empNorm = xc_empNorm ./ stdPerm;
% 
% %% VISUALIZE XC
% figure(2); clf; hold on;
% plot(lags, xc_empNorm, 'linew', 3, 'color', [0 0 0])
% ax            = gca;
% ax.FontSize   = 20;
% ax.FontWeight = 'bold';
% box off
% title('Crosscorrelation between encoding and retrieval spiking (0s : 3s)')
% xlabel('Lags [ms] | Positive peak indicates earlier retrieval')
% ylabel('Z-Score (normalized with permutation)')
% plot([0 0], get(gca, 'Ylim'), '--', 'color', [0 0 0 0.5], 'linew', 2)
%     
% figure;
% plot(lags, xc_emp, 'Color', [0 0 1]);hold on
% plot(lags, meanPerm, 'Color', [1 0 0]);legend Emp Rand;hold off
% 
% figure;
% hist(rlocs,20);
% hold on
% plot([eloc eloc],[0 100],'--', 'Color', [1 0.2 0.2], 'Linewidth',1);
% pval=numel(find(rlocs>=eloc))/nperm;
% title(['Pval = ', num2str(pval)]);

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