% %% PERMUTATION chanMean
% % clear all
% load('X:\Luca\data\allSbj\preCuePowDiff_orthNorm.mat', 'bundlePow_chanMean', 'hz')
% 
% % bundlePow_chanMean.idx is already the channel mean in the form of
% % trl x power
% % so here I average over trials and 
% % then concatenate all bundles
% 
% % TAKE THE MEDIAN OVER ALL TRIALS FOR EACH BUNDLE ("{bundles}: 1 x pow")
% idxPow = cellfun(@(x) median(x,1), bundlePow_chanMean.idx, 'un', 0 ); 
% idxPow = vertcat(idxPow{:}); % "bundles x pow"
% 
% % TAKE THE MEDIAN OVER ALL BUNDLES 
% idxPow = nanmedian(idxPow,1);
% 
% % same for ndx
% ndxPow = cellfun(@(x) median(x,1), bundlePow_chanMean.ndx, 'un', 0 );
% ndxPow = vertcat(ndxPow{:});
% ndxPow = nanmedian(ndxPow,1);
% 
% % DIFFERENCE BETWEEN BUNDLES
% difPow = idxPow - ndxPow;
% 
% %% Note: median(a) - median(b) ~= median(a-b)
% %  Also: In a AxBxC matrix it matters which dimension you median first
% %  (although not much)
% 
% nperm = 1000;
% difPowPerm = zeros(nperm, size(difPow,2));
% for perm = 1 : nperm
%     
%     difPowPermBund = [];
%     for bund = 1 : size(idxPow,1) % LOOP OVER EACH BUNDLE
%         numIdx = size(bundlePow_chanMean.idx{bund}, 1);                        % number of indexed trials in that bundle
%         allPow = [bundlePow_chanMean.idx{bund}; bundlePow_chanMean.ndx{bund}]; % "trial x pow" for idx and ndx trials concatenated
%         allPow = allPow(randperm (size(allPow,1))  ,: );                       % shuffle "trial x pow" matrix
%         
%         idxPowPerm = median(allPow(1:numIdx,:),1); 
%         ndxPowPerm = median(allPow(numIdx+1:end,:),1);
%         
%         difPowPermBund(bund,:) = idxPowPerm - ndxPowPerm;
%         
%     end
%     difPowPerm(perm,:) = nanmedian(difPowPermBund,1);
% end
% 
% 
% 
% figure(1); clf; hold on;
% plot(difPow, 'linew',2);
% hiTH = prctile(difPowPerm,100-2.5,1);
% loTH = prctile(difPowPerm,0+2.5,1);
% plot(hiTH, 'color', 'r');
% plot(loTH, 'color', 'r');
% xlim([1 100])
% ylim([-0.01 0.001])
% 
% 
% 
% %% FREQUENCY BANDS
% freq.delta  = hz<3;
% freq.theta  = hz>= 3 & hz< 8;
% freq.alphaBeta = hz>= 8 & hz<20;
% freq.gamma = hz>=40 & hz<100; 
% ntest = 4;
% 
% %% BINNED
% difPowBin     = [mean(difPow(:,freq.delta),2) mean(difPow(:,freq.theta),2) mean(difPow(:,freq.alphaBeta),2) mean(difPow(:,freq.gamma),2)];   
% difPowPermBin = [mean(difPowPerm(:,freq.delta),2) mean(difPowPerm(:,freq.theta),2) mean(difPowPerm(:,freq.alphaBeta),2) mean(difPowPerm(:,freq.gamma),2)];   
% 
% hiTH = prctile(difPowPermBin,100-(2.5/ntest),1);
% loTH = prctile(difPowPermBin,0+(2.5/ntest),1);
% 
% % Visualize FREQUENCYBANDS
% plot([1    3], [-0.009  -0.009 ], 'linew', 2, 'color', 'k');
% plot([3    8], [-0.00925 -0.00925], 'linew', 2, 'color', 'k');
% plot([8   20], [-0.009  -0.009 ], 'linew', 2, 'color', 'k');
% plot([40 100], [-0.00925 -0.00925], 'linew', 2, 'color', 'k');
% 
% 
% xticks([2 5.5 14 70])
% xticklabels({'Delta', 'Theta', 'Alpha Beta', 'Gamma'})
% ylabel('Power Indexed Trials Minus Non-Indexed Trials')
% yticks('')
% 
% scatter(2,   -0.0085, 100, 'k', 'filled', 'p')
% scatter(5.5, -0.0085, 100, 'k', 'filled', 'p')
% scatter(14,  -0.0085, 100, 'k', 'filled', 'p')
% scatter(70,  -0.0085, 100, 'k', 'filled', 'p')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% QUESTIONS:
%  median over channel because of noisy channels?
%  don't mean over channel because they are independet?

%  mean [idx(chan x pow) - ndx(chan x pow) ]
%  instead of
%  mean[idx(chan x pow)] - mean [ndx(chan x pow)]
% ==> mean the channel difference instead of meaning the idx/ndx power and
% then substract

%  mean instead median over bundles because they are independent??

% Für die Permutation jedoch müsste ich nun innerhalb eines jeden Channels die
% trial definition (idx / ndx) shufflen und anschließend erst das channelspezifische
% Powerspektrum berechnen und anschließend über alle channel averagen.



%% EMPIRICAL DIFFERENCE
clear
load('X:\Luca\data\allSbj\preCuePowDiff_orthNorm.mat', 'bundlePow_allChan', 'hz');
nperm = 1000;
difPerm = zeros(nperm, 797);
%  MEDIAN OVER TRIAL
%  SUBSTRACT NDX FROM IDX
%  MEAN OVER CHANNEL
%  MEAN OVER BUNDLES
idxPow = cellfun(@(x) median(x, 2), bundlePow_allChan.idx, 'un', 0); % MEDIAN OVER TRIALS
ndxPow = cellfun(@(x) median(x, 2), bundlePow_allChan.ndx, 'un', 0); % MEDIAN OVER TRIALS

diffPow = cellfun(@minus, idxPow, ndxPow, 'un', 0);     % SUBSTRACT NDX FROM IDX
diffPow = cellfun(@(x) mean(x, 1), diffPow, 'un', 0);   % MEDIAN OVER CHANNEL (" 1 x 1 x POW ")
diffPow = cellfun(@squeeze, diffPow, 'un', 0);          % SQUEEZE DIMENSION
diffPow = cat(2, diffPow{:})';                          % CONCATENATE OVER BUNDLES

%% Visu
figure(1); clf;
subplot(211); hold on;
title('Mean over Channel')
for ii = 1:size(diffPow,1); plot(diffPow(ii,:), 'color', [0 0 0 0.05], 'linew', 5);end
xlim([1 75])
subplot(212); hold on;
title('Median or Mean over Bundles')
plot(nanmedian(diffPow,1), 'linew', 2, 'color', 'b');
plot(nanmean(diffPow,1),   'linew', 2, 'color', 'k');
xlim([1 75])
ylim([-0.0020 0.0160])

diffPow = nanmean(diffPow, 1);                      % MEDIAN OVER BUNDLES

difPerm = zeros(nperm, 797);
for perm = 1 : nperm
    tic
    disp(perm)
    
    shufDiffBund = zeros(94,797);
    for bund = 1:94
        numIdx   = size(bundlePow_allChan.idx{bund},2);                              % number of indexed trials in that bundle (" chan x trl x pow ")
        shufPow  = cat(2, bundlePow_allChan.idx{bund}, bundlePow_allChan.ndx{bund}); % concatenate idx&ndx trials
        randIdx  = randperm(size(shufPow, 2)); 
        shufPow  = shufPow(:,randIdx,:);                                             % shuffle "chan x trl x pow" matrix
        
        shufIdx  = shufPow(:,1:numIdx,:); % EXTRACT "CHAN x TRL x POW" FOR idx
        shufIdx  = squeeze(median(shufIdx,2)); % MEDIAN OVER TRIALS
        
        shufNdx  = shufPow(:,numIdx+1:end,:);
        shufNdx  = squeeze(median(shufNdx,2)); % MEDIAN OVER TRIALS
        
        shufDiff = shufIdx - shufNdx; 
        
        shufDiffBund(bund,:) = nanmean(shufDiff,1);
    end
difPerm(perm,:) = nanmean(shufDiffBund,1); 
toc
end

%% Visu
figure(2); clf; hold on;
title('Difference PreCue Power')
plot(hz, nanmean(diffPow,1),   'linew', 2, 'color', 'b');
plot(hz, prctile(difPerm, 95, 1), 'linew', 1.5, 'color', 'r');
plot(hz, prctile(difPerm,  5, 1), 'linew', 1.5, 'color', 'r');
xlim([1 75])
xlabel('Freqzency [Hz]')
ylabel('Power (Idx-Ndx)')
% ylim([-0.0020 0.0160])

