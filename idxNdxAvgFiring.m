% average rate: index units
% average rate: ndx units

clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')

encRet = 'enc';
% dt = -1:0.001:10;
dt = -1.375:0.001:10.375+0.001;
dt = dt-0.0005; % center around 0


%% gaussian kernel
mlength = [-0.075:0.0002:0.075];
mSigma  = 0.02; % 20ms
% mSigma  = 0.01;
mKernel = normpdf(mlength,0,mSigma);
mKernel = mKernel/max(mKernel); % normalize peak to 1


ndx  = [];
idxN = [];
idxI = [];
for su = 1: length(allSpks)
    
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
    
    
    %% NORMALISATION OVER WHOLE TRIAL
%     convBL   = mean(allConv,2);
%     convMean = mean(convBL);
%     convSTD  = std(convBL,0,1);
%     
%     allConv  = allConv-convMean;
%     allConv  = allConv./convSTD;
    
    
%     %% ALTERNATIVE NORMALISATION OVER PRECUE PERIOD
    convBL   = mean(allConv(:,1:1001),2);
    convMean = mean(convBL);
    convSTD  = std(convBL,0,1);
    
    allConv  = allConv-convMean;
    allConv  = allConv./convSTD;
    
    switch sum(idxTrl)>0
        case 0 % ndx unit
            ndx = [ndx; mean(allConv,1)];
        case 1 % idx unit
            idxN = [idxN; mean(allConv(~idxTrl,:),1)];
            idxI = [idxI; mean(allConv( idxTrl,:),1)];
    end
    fprintf(repmat('\b',1,length(lineLength)-1));

end



%% VISUALISATION
figure(1); clf;
subplot(1,2,1); hold on;
plotDT = [-1000:1:10000];
plotDat = nanmedian(ndx,1); %% SOME DO NOT FIRE DURING PRECUE PERIOD!!
plot(plotDT, plotDat, 'color', [0 0 0], 'linew', 2)

shade     = std(ndx,0,1) / sqrt(size(ndx,1));
shadeDT   = [plotDT, fliplr(plotDT)];
inBetween = [plotDat-shade, fliplr(plotDat+shade)];
fillHand = fill(shadeDT, inBetween, [0.1059, 0.6196, 0.4667] );
fillHand.FaceAlpha = 0.5;


xlim([-1000 10000])
ylim([-0.5 3.5])
xlabel('Time')
xticks([0 2000 5000 10000])
switch strcmp(encRet, 'enc')
    case 0
        xticklabels({'Cue', 'Q', '5s','10s'})
    case 1
        xticklabels({'Cue', 'Stim', '5s','10s'})
end
yticks(0)
ylabel('Firing rate (z-values)')
title('Non Index Units')
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
box off

%%
% IU INDEXED TRIALS
subplot(1,2,2); hold on;
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
yticks(0)
ylabel('Firing rate (z-values)')
title('Index Units')

% % LEGENDS
    L(1) = plot(nan, nan, 'color', [0.4588, 0.4392, 0.7020]);
    L(2) = plot(nan, nan, 'color', [0.1059, 0.6196, 0.4667]);
[legPos, hobj, ~, ~] = legend(L, {'Indexed Trials', 'Non-indexed Trials'}, 'FontSize',16, 'FontWeight', 'bold');
set(hobj,'LineWidth',15);
set(legPos, 'Position', [0.7629    0.7858    0.1962    0.1022])
legend('boxoff')

% hobj = legend('a')
% c=get(hobj,'Children')
% set(c,'Color',[0 0 1])

%     legend('boxoff');

switch strcmp(encRet, 'enc')
    case 0
        sgtitle('Retrieval', 'FontSize', 20, 'FontWeight', 'bold')
    case 1
        sgtitle('Encoding',  'FontSize', 20, 'FontWeight', 'bold')
end

ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
box off



     
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
