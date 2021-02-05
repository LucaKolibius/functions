%% angles:
% play around with TOI
% play around with kernel sigma
% play around with mean within SU

function [allConvEnc, allConvRet] = xcEncRet
whereAmI(0);
global prePath;
load([prePath, 'Luca\data\allSbj\allSpksHZ.mat'], 'allSpks')

allConvEnc = [];
allConvRet = [];

%% make spike density
for su = 1:size(allSpks,2)
    disp(su);
    
    if sum(allSpks(su).idxTrl) == 0
        continue % NO INDEX TRIALS
    end
    
    [hitsIdx, ~, spksEnc, spksRet, ~, ~, ~, ~, ~, ~] = loadInDat(allSpks, su);
    spksEnc = spksEnc'; spksRet = spksRet';
    
  
    %% Visualization
    timeWindow = [-0.375 2.375];
    dt         = linspace(timeWindow(1), timeWindow(2), (abs(timeWindow(1)) + abs(timeWindow(2))) *1000+1 );
    ploDT     = dt+0.0005; ploDT(end) = [];

    
    % % Encoding
    n_encHit = [];
    counter1 = 0;
    for trl = 1:size(hitsIdx,1) % trials
        counter1 = counter1+1;
        x = spksEnc{1, hitsIdx(trl)};
        [n_encHit(counter1,:),~] = histcounts(x,dt);
    end
    
 
    % % RETRIEVAL
    counter1 = 0;
    n_retHit = [];
    for trl = 1 : size(hitsIdx,1) % trials
        counter1 = counter1+1;
        x = spksRet{1, hitsIdx(trl)};
        [n_retHit(counter1,:),~] = histcounts(x,dt);
    end
    
    
    retIdxSpikes = n_retHit(allSpks(su).idxTrl,:);
    retIdxSpikes = sum(retIdxSpikes,1);
    encIdxSpikes = n_encHit(allSpks(su).idxTrl,:);
    encIdxSpikes = sum(encIdxSpikes,1);
    
    % gaussian kernel
    mlength = [-0.075:0.0002:0.075];
    mSigma  = 0.02; % 20ms
    % mSigma  = 0.01;
    mKernel = normpdf(mlength,0,mSigma);
    mKernel = mKernel/max(mKernel); % normalize peak to 1
    
    % encoding
    encConv   = conv(mKernel, encIdxSpikes);
    % get rid of edges
    encConv(1:750)       = [];
    encConv(end-749:end) = [];
    
    
        
    % retrieval
    retConv = conv(mKernel,retIdxSpikes);
    % get rid of edges
    retConv(1:750)       = [];
    retConv(end-749:end) = [];
    
    allConvEnc = [allConvEnc; encConv];
    allConvRet = [allConvRet; retConv];
    
%     % plot convolved series
%     hold on
%     title('Firing rate during indexed trials')
%     plot(ploDT, encConv, 'linew', 3, 'color', mBlue);
%     plot(ploDT, retConv, 'linew', 3, 'color', mOrange);
%     ylabel({'Smoothed', 'Occurences', ''})
%     box on
%     xlabel('Time after cue onset (seconds)')
%     xticks([-2 0 2 5:5:20]);
%     yticks('')
%     xlim([-2 20])
% 
%     xticklabels({'-2', 'Cue Onset', '2', '5', '10', '15', '20'})
%     
%     mAx = gca;
%     mAx.YAxis.FontWeight = 'bold';
%     mAx.XAxis.FontWeight = 'bold';
%     mAx.FontSize = 14;
%     line([0 0],[get(gca, 'YLim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3) % cue onset
%     line([2 2],[get(gca, 'YLim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3) % stimulus onset
%     
%     hold on
%     L(1) = plot(nan, nan, 'color', mBlue);
%     L(2) = plot(nan, nan, 'color', mOrange);
%     [legPos, hobj, ~, ~] = legend(L, {'Encoding', 'Retrieval'}, 'FontSize',16, 'FontWeight', 'bold');
%     hl = findobj(hobj,'type','line');
%     set(hl,'LineWidth',15);
%     set(legPos, 'Position', [0.033 0.875 0.078 0.061])
%     legend('boxoff');


end % end of loop


%% XC
XC = [];
locsEnc = [];
locsRet = [];
for itrl = 1:size(allConvEnc,1)
    maxlag = 1000;
%     [XC(itrl,:),LAGS] = xcorr(allConvEnc(itrl,:), allConvRet(itrl,:), maxlag);
    [XC(itrl,:),LAGS] = crosscorr(allConvEnc(itrl,:), allConvRet(itrl,:), 'NumLags', maxlag);

    [~,allpks] = findpeaks(allConvEnc(itrl,:));
    if ~isempty(allpks)
        locsEnc(itrl) = allpks(1);
    else
        locsEnc(itrl) = NaN;
    end
    
    [~,allpks] = findpeaks(allConvRet(itrl,:));
    if ~isempty(allpks)
        locsRet(itrl) = allpks(1);
    else
        locsRet(itrl) = NaN;
    end
    
end

%% location of first peak for Enc and Ret
nanidx = isnan(locsRet)| isnan(locsEnc);
locsRet(nanidx) = NaN;
locsEnc(nanidx) = NaN;

figure(3); clf;
dt = 0:100:2000;
histogram(locsEnc,dt); hold on;
histogram(locsRet,dt);

%% average crosscorrelation (meaned within IU)
figure(1); clf;hold on;
plot(LAGS, nanmean(XC,1), 'linew', 3);
plot([0 0],get(gca, 'Ylim'), '--', 'color', [0 0 0 0.5], 'linew', 2)
xlabel('Lag')
ylabel('r')
mAx = gca;
mAx.FontWeight = 'bold';
mAx.FontSize = 20;

%% grand mean spike density for enc and ret
figure(2); clf;
plot(mean(allConvEnc,1), 'color', [0.368627450980392,0.235294117647059,0.600000000000000], 'linew', 2); hold on;
plot(mean(allConvRet,1), 'color', [0.901960784313726,0.380392156862745,0.00392156862745098], 'linew', 2)
mAx = gca;
xlabel('Time [ms]')
ylabel('Firing Density')
legend({'Encoding','Retrieval'})
mAx.FontWeight = 'bold';
mAx.FontSize = 20;


%% crosscorrelation of grand mean spike density enc X ret
[xcf,lags,bounds]=crosscorr(mean(allConvEnc,1), mean(allConvRet,1), 'NumLags', maxlag);
plot(lags, xcf, 'linew', 3)

%% find peak distance between grand averages
enc = mean(allConvEnc,1);
ret = mean(allConvRet,1);

[~,locEnc] = findpeaks(enc);
locEnc     = locEnc(1);
[~,locRet] = findpeaks(ret);
locRet     = locRet(1);
peakDist   = locRet-locEnc; %% <-- empirical distance between first peaks

%% shuffled distribution of peak distance between enc and ret grand averages
allConv = [allConvEnc; allConvRet];
nperm = 1000;
peakDistPerm = [];
for perm = 1 : nperm
    disp(perm)
    allConv = allConv(randperm(size(allConv,1)),:);
    permEnc = allConv(1:size(allConv,1)/2,:);
    permRet = allConv(size(allConv,1)/2+1:end,:);
    
    permEnc = mean(permEnc,1);
    permRet = mean(permRet,1);
    
    [~,locEnc]         = findpeaks(permEnc);
    locEnc             = locEnc(1);
    [~,locRet]         = findpeaks(permRet);
    locRet             = locRet(1);
    peakDistPerm(perm) = locRet-locEnc;
end

prctile(peakDistPerm,95)
prctile(peakDistPerm,5)
histogram(peakDistPerm) %% <-- seems like there are two clusters 


end % end of function