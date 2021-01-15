function preCuePowDiff_perm(highLow)

switch highLow
    case 1 % HIGH
        load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\preCuePowDiff_orthDeMea.mat', 'allSUPowHigh');
        allSUPow = allSUPowHigh;
        hz = 70:1:150;
    case 2 % LOW
        load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\preCuePowDiff_orthDeMea.mat', 'allSUPowLow');
        allSUPow = allSUPowLow;
        hz = 1:1:45;     
end

%  MEDIAN OVER TRIAL
%  SUBSTRACT NDX FROM IDX
%  MEAN OVER BUNDLES
idxPow = cellfun(@(x) nanmedian(x, 1), allSUPow.idx, 'un', 0); % MEDIAN OVER TRIALS
ndxPow = cellfun(@(x) nanmedian(x, 1), allSUPow.ndx, 'un', 0); % MEDIAN OVER TRIALS

diffPow = cellfun(@minus, idxPow, ndxPow, 'un', 0);         % SUBSTRACT NDX FROM IDX
% diffPow = cellfun(@(x) mean(x, 1), diffPow, 'un', 0);       % MEDIAN OVER CHANNEL (" 1 x 1 x POW ")
% diffPow = cellfun(@squeeze, diffPow, 'un', 0);              % SQUEEZE DIMENSION
diffPow = cat(1, diffPow{:});                              % CONCATENATE OVER SU

%% Visu
figure(1); clf;
subplot(211); hold on;
title('Mean over Channel')
for ii = 1:size(diffPow,1); plot(hz, diffPow(ii,:), 'color', [0 0 0 0.05], 'linew', 5);end
% xlim([1 80])
subplot(212); hold on;
title('Median or Mean over Bundles')
plot(hz, nanmedian(diffPow,1), 'linew', 2, 'color', 'b');
plot(hz, nanmean(diffPow,1),   'linew', 2, 'color', 'k');
% xlim([1 80])
% ylim([-0.0020 0.0160])
legend({'median', 'mean'})
% plot([80 140],[0 0 ], 'k--');

diffPow = nanmean(diffPow, 1);                           % MEDIAN OVER SU

nperm = 10000;
numSU = size(allSUPow.idx,2);
difPerm = zeros(nperm, length(hz));

for perm = 1 : nperm
    tic
    disp(perm)
    
        shufDiffSU = zeros(numSU,length(hz));
    for su = 1:numSU
        numIdx   = size(allSUPow.idx{su},1);                   % number of trials that this SU indexed
        shufPow  = cat(1, allSUPow.idx{su}, allSUPow.ndx{su}); % concatenate idx&ndx trials
        randIdx  = randperm(size(shufPow, 1)); 
        shufPow  = shufPow(randIdx, :);                        % shuffle "trl x pow" matrix
        
        % IDX
        shufIdx  = shufPow(1:numIdx,:);                        % EXTRACT "TRL x POW" FOR idx
        shufIdx  = nanmedian(shufIdx,1);                          % MEDIAN OVER TRIALS
        
        % NDX                                                  & REPEAT FOR NDX
        shufNdx  = shufPow(numIdx+1:end,:);
        shufNdx  = nanmedian(shufNdx,1);                 
        
        shufDiff = shufIdx - shufNdx; 
        
        shufDiffSU(su,:) = shufDiff;
    end
difPerm(perm,:) = nanmean(shufDiffSU,1);                        % MEDIAN OVER SU
toc
end

%% Visu (NON-BINNED)
figure(2); clf; hold on;
title('Difference Ripple Power -  Whole Trial (85-140 hz)')
plot(hz, diffPow,   'linew', 2, 'color', 'b');
plot(hz, prctile(difPerm, 100-5/6, 1), 'linew', 1.5, 'color', 'r');
plot(hz, prctile(difPerm,  5/6, 1), 'linew', 1.5, 'color', 'r');
xlim([hz(1) hz(end)])
xlabel('Freqzency [Hz]')
ylabel('Power (Idx-Ndx)')
legend({'Power Difference', 'Upper & Lower Threshold (corrected)'});
% ylim([-0.0020 0.0160])


%% Visu (BINNED)
switch highLow
    case 1 % HIGH
      
        freq.lowRip  = hz >= 80  & hz <= 100;
        freq.medRip  = hz >= 100 & hz <= 120;
        freq.highRip = hz >= 120 & hz <= 140;
        
        
        
        powDiffBin   = [sum(diffPow(:,freq.lowRip),2) sum(diffPow(:,freq.medRip),2) sum(diffPow(:,freq.highRip),2)];
        diffPermBin  = [sum(difPerm(:,freq.lowRip),2) sum(difPerm(:,freq.medRip),2) sum(difPerm(:,freq.highRip),2)];
        
        figure(3); clf; hold on;
        title('Difference PreCue Power')
        plot(powDiffBin, 'linew', 2, 'color', 'b');
        plot(prctile(diffPermBin, 100-5/3, 1), 'linew', 1.5, 'color', 'r');
        plot(prctile(diffPermBin,  5/3, 1), 'linew', 1.5, 'color', 'r');
        xlabel('Frequency [Hz]')
        ylabel('Power (Idx-Ndx)')
        xticks([1:3])
        xticklabels({'Slow Ripple (80-100 hz)', 'Medium Ripple (100-120 hz)', 'Fast Ripple (120-140 hz)'})
        
    case 2 % LOW
        freq.delta  = hz >= 2.5 & hz <=  5;
        freq.theta  = hz >= 5.5 & hz <= 10;
        freq.alphaL = hz >= 10  & hz <= 19;
        freq.alphaH = hz >= 19  & hz <= 32;
        freq.beta   = hz >= 34  & hz <= 100;
        
        powDiffBin  = [sum(diffPow(:,freq.delta),2) sum(diffPow(:,freq.theta),2) sum(diffPow(:,freq.alphaL),2) sum(diffPow(:,freq.alphaH),2) sum(diffPow(:,freq.beta),2)];
        diffPermBin = [sum(difPerm(:,freq.delta),2) sum(difPerm(:,freq.theta),2) sum(difPerm(:,freq.alphaL),2) sum(difPerm(:,freq.alphaH),2) sum(difPerm(:,freq.beta),2)];
        
        figure(3); clf; hold on;
        title('Difference PreCue Power')
        plot(powDiffBin,   'linew', 2, 'color', 'b');
        plot(prctile(diffPermBin, 100-5/5, 1), 'linew', 1.5, 'color', 'r');
        plot(prctile(diffPermBin,  5/5, 1), 'linew', 1.5, 'color', 'r');
        xlabel('Frequency [Hz]')
        ylabel('Power (Idx-Ndx)')
        xticks([1:5])
        xticklabels({'Slow Theta (2.5-5)', 'Fast Theta (5.5-10)', 'Alpha (10-19)', 'Beta (19-32)', 'Gamma (34-100)'})
        % ylim([-0.0020 0.0160])
        
end


%% POOLED OVER IDX & NDX
idx = vertcat(allSUPow.idx{:});
ndx = vertcat(allSUPow.ndx{:});

figure(1); clf; 
subplot(211); hold on;
title('Idx')
for ii = 1:size(idx,1); plot(hz, idx(ii,:), 'color', [0 0 0 0.05], 'linew', 5);end
% xlim([1 80])
subplot(212); hold on;
title('Ndx')
for ii = 1:size(ndx,1); plot(hz, ndx(ii,:), 'color', [0 0 0 0.05], 'linew', 5);end

diffPowMean = nanmean(idx,1) - nanmean(ndx,1);
diffPowMedian = nanmedian(idx,1) - nanmedian(ndx,1);
figure(2); clf; hold on;
plot(hz, diffPowMean, 'color', 'b');
plot(hz, diffPowMedian, 'color', 'r');
diffPow = diffPowMean;

figure
plot(hz, nanmedian(diffPow,1), 'linew', 2, 'color', 'b');
plot(hz, nanmean(diffPow,1),   'linew', 2, 'color', 'k');

nperm = 10000;
for perm = 1 : nperm
    tic
    disp(perm)
    
    numIdx   = size(idx,1);                   % number of trials that this SU indexed
    shufPow  = [idx; ndx]; % concatenate idx&ndx trials
    randIdx  = randperm(size(shufPow, 1));
    shufPow  = shufPow(randIdx, :);                        % shuffle "trl x pow" matrix
    
    % IDX
    shufIdx  = shufPow(1:numIdx,:);                        % EXTRACT "TRL x POW" FOR idx
    
    % NDX                                                  & REPEAT FOR NDX
    shufNdx  = shufPow(numIdx+1:end,:);
    
    shufDiff = nanmean(shufIdx) - nanmean(shufNdx);
    
    diffPerm(perm,:) = shufDiff;
    toc
end

figure(3); clf; hold on;
plot(hz, diffPow, 'color', 'k', 'linew', 2);
plot(hz, prctile(diffPerm, 100-5, 1), 'linew', 1.5, 'color', 'r');
plot(hz, prctile(diffPerm,  5, 1), 'linew', 1.5, 'color', 'r');


end % OF FUNCTION
