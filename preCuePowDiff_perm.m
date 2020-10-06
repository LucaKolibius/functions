%% PERMUTATION chanMean
% clear all
load('X:\Luca\data\allSbj\preCuePowDiff_orthNorm.mat', 'bundlePow_chanMean', 'hz')

% bundlePow_chanMean.idx is already the channel mean in the form of
% trl x power
% so here I average over trials and 
% then concatenate all bundles
idxPow = cellfun(@(x) mean(x,1), bundlePow_chanMean.idx, 'un', 0 );
idxPow = vertcat(idxPow{:}); 
% idxPow = nanmean(idxPow,1);   % MEAN OF INDEX
idxPow = nanmedian(idxPow,1); % MEDIAN OF INDEX

% same for ndx
ndxPow = cellfun(@(x) mean(x,1), bundlePow_chanMean.ndx, 'un', 0 );
ndxPow = vertcat(ndxPow{:});
% ndxPow = nanmean(ndxPow,1);   % MEAN OF NDX
ndxPow = nanmedian(ndxPow,1); % MEDIAN OF NDX

% difference between bundles
difPow = idxPow - ndxPow;

nperm = 1000;
difPowPerm = zeros(nperm, size(difPow,2));
for perm = 1 : nperm
    
    difPowPermBund = [];
    for bund = 1 : size(idxPow,1)
        numIdx = size(bundlePow_chanMean.idx{bund}, 1);
        allPow = [bundlePow_chanMean.idx{bund}; bundlePow_chanMean.ndx{bund}];
        allPow = allPow(randperm (size(allPow,1))  ,: );
        
        idxPowPerm = mean(allPow(1:numIdx,:),1);
        ndxPowPerm = mean(allPow(numIdx+1:end,:),1);
        
        difPowPermBund(bund,:) = idxPowPerm - ndxPowPerm;
        
    end
    difPowPerm(perm,:) = nanmedian(difPowPermBund,1);
end



figure(1); clf; hold on;
plot(difPow, 'linew',2);
hiTH = prctile(difPowPerm,100-2.5,1);
loTH = prctile(difPowPerm,0+2.5,1);
plot(hiTH, 'color', 'r');
plot(loTH, 'color', 'r');
xlim([1 100])
ylim([-0.01 0.001])



%% FREQUENCY BANDS
freq.delta  = hz<3;
freq.theta  = hz>= 3 & hz< 8;
freq.alphaBeta = hz>= 8 & hz<20;
freq.gamma = hz>=40 & hz<100; 
ntest = 4;

%% BINNED
difPowBin     = [mean(difPow(:,freq.delta),2) mean(difPow(:,freq.theta),2) mean(difPow(:,freq.alphaBeta),2) mean(difPow(:,freq.gamma),2)];   
difPowPermBin = [mean(difPowPerm(:,freq.delta),2) mean(difPowPerm(:,freq.theta),2) mean(difPowPerm(:,freq.alphaBeta),2) mean(difPowPerm(:,freq.gamma),2)];   

hiTH = prctile(difPowPermBin,100-(2.5/ntest),1);
loTH = prctile(difPowPermBin,0+(2.5/ntest),1);

% Visualize FREQUENCYBANDS
plot([1    3], [-0.009  -0.009 ], 'linew', 2, 'color', 'k');
plot([3    8], [-0.00925 -0.00925], 'linew', 2, 'color', 'k');
plot([8   20], [-0.009  -0.009 ], 'linew', 2, 'color', 'k');
plot([40 100], [-0.00925 -0.00925], 'linew', 2, 'color', 'k');


xticks([2 5.5 14 70])
xticklabels({'Delta', 'Theta', 'Alpha Beta', 'Gamma'})
ylabel('Power Indexed Trials Minus Non-Indexed Trials')
yticks('')

scatter(2,   -0.0085, 100, 'k', 'filled', 'p')
scatter(5.5, -0.0085, 100, 'k', 'filled', 'p')
scatter(14,  -0.0085, 100, 'k', 'filled', 'p')
scatter(70,  -0.0085, 100, 'k', 'filled', 'p')