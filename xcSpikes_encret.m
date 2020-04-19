%% xc between cue and cue+1 spiketimes enc+ret
%load in data
cd('X:\Luca\indexSUvisu\varWindow_ENC_-1-113_RET_-1113_DP_p99_th1')
clear
load('globalTest')
nperm = 10000;

%%
% gaussian kernel
mlength = [-0.075:0.0002:0.075];
mSigma  = 0.02; % 20ms
% mSigma  = 0.01;
mKernel = normpdf(mlength,0,mSigma);
mKernel = mKernel/max(mKernel); % normalize peak to 1
dt = [0:1:1000];

% if there are no spikes either during encoding or retrieval, kick trial out
insuffSpks = or(cellfun(@isempty, encVec_allSU), cellfun(@isempty, retVec_allSU)); % 1 if no spikes between cue and cue+1
encVec_allSU(insuffSpks) = [];
retVec_allSU(insuffSpks) = [];

allEnc = [];
allRet = [];
for iEnc = 1:size(encVec_allSU,2)
    trlN = encVec_allSU{iEnc}*1000;  % spiketimes for first encoding trial in ms
    binVal = histcounts(trlN, dt);   % move into 1ms bins ranging from 0 to 1000ms
    encConv = conv(mKernel, binVal); % convolute with a gaussian
    
    encConv(1:375)       = [];  % cut off wings
    encConv(end-374:end) = [];  % cut off wings
    allEnc = [allEnc, {encConv}]; % save that spiketime series
end

% same for ret
for iRet = 1:size(retVec_allSU,2)
    trlN = retVec_allSU{iRet}*1000;
    binVal = histcounts(trlN, dt);
    retConv = conv(mKernel, binVal);
    
    retConv(1:375)       = [];
    retConv(end-374:end) = [];
    allRet = [allRet, {retConv}];
end

%% visualize original series
close all
    mBlue = [0 0.4470 0.7410];
    mOrange = [0.8500 0.3250 0.0980];

figure
hold on
plot([allEnc{:}], 'linewidth', 3); plot([allRet{:}], 'linewidth', 3); % plot spiketime series
% legend
    hold on
    L(1) = plot(nan, nan, 'color', mBlue);
    L(2) = plot(nan, nan, 'color', mOrange);
    [legPos, hobj, ~, ~] = legend(L, {'Encoding', 'Retrieval'}, 'FontSize',18, 'FontWeight', 'bold');
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',3);
    set(legPos, 'Position', [0.79,0.837,0.11,0.09])
    legend('boxoff');
    
    xticks(500:1000:9000)
    xticklabels(string(1:9))
    xlabel('Indexed Trial Number')
    yticks('')
    ylabel('Smoothed Firing')
    
    
    for j = 1:9
    line([1000*j 1000*j],[get(gca,'Ylim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',2)
    end
    
    mAx = gca;
    mAx.FontWeight = 'bold';
    mAx.FontSize = 18;

% compute empirical crosscorrelation
[a,b] = crosscorr([allEnc{:}],[allRet{:}], 1000); % crosscorrelate encoding and retrieval spiketime series

%% permutation
numIdxTrl = size(retVec_allSU,2); % number of index trials that we are looking into
permRet = zeros(nperm, 1000*numIdxTrl ); % preallocation
lags = 1000;
permXC = zeros(nperm, lags*2+1); % preallocation

for perm = 1:nperm
%     disp(perm)
    randIdx = randperm(numIdxTrl); % random index
    permRet(perm,:)  = [allRet{randIdx}]; % scramble retrieval trials
    [permXC(perm,:), lagVec] = crosscorr([allEnc{:}], permRet(perm,:), lags);
end
%% in this version allEnc and allRet were 1x3000 doubles instead of 3 cells with 1x1000 doubles in them (as is now)
% numIdxTrl = size(retVec_allSU,2); % number of index trials that we are looking into
% permRet = zeros(nperm, 1000*numIdxTrl ); % preallocation
% for perm = 1:nperm
%     randIdx = randperm(numIdxTrl); % random index
%     newRet  = retVec_allSU(randIdx); % scramble retrieval trials
%     
%     temp = [];
%     for ib = 1:numIdxTrl
%         trlN = newRet{ib}*1000;
%         binVal = histcounts(trlN, dt);
%         
%         permConv = conv(mKernel, binVal);
%         permConv(1:375)       = [];
%         permConv(end-374:end) = [];
%         
%         temp = [temp, permConv];
%     end
%     
%     permRet = [permRet; temp];
% end
% 
% lags = 1000;
% permXC = zeros(nperm, lags*2+1);
% for xc = 1:nperm
%     [permXC(xc,:), lagVec] = crosscorr(allEnc, permRet(xc,:), lags);
% end

%% 
% avgXC = mean(permXC,1);
% stdXC = std(permXC,1);
top = prctile(permXC, 95);
bot = prctile(permXC, 5);

% plot xcor
figure
hold on
plot(b,a, 'linewidth', 3 ) % max at 194

xlabel('Lag between Encoding and Retrieval (ms)');
title('Earlier firing during retrieval ({\it r} = 0.58 @ 194ms )')
ylabel('Correlation')

mAx = gca;
mAx.FontWeight = 'bold';
mAx.FontSize = 18;
ylim([-1 1])
plot_dt = -1000:1:1000;
plot(plot_dt, top, '--', 'linewidth', 3, 'color', 'r' )
plot(plot_dt, bot, '--', 'linewidth', 3, 'color', 'r' )
line([0 0],[get(gca,'Ylim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3)
