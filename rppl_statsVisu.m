%% ANALYSIS PART
clear
% load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks'); % I no longer do spike ripple coincidence
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_rppls.mat', 'allSpks');
% numMW = 1; % for macro recording
tw = 0; % :1:3 % only whole trial

        numRpplIdx = [];
        numRpplNdx = [];
        lenRpplNdx = [];
        lenRpplIdx = [];
        rpplDenIdx = [];
        rpplDenNdx = [];
%         coincIdx   = [];
%         coincNdx   = [];

for su = 1 : length(allSpks)
    
    disp(su)
    if any(isnan(allSpks(su).idxTrlSing))
        continue
    end
    
    bidsID     = allSpks(su).bidsID;
    sesh       = allSpks(su).sesh;
    bundlename = allSpks(su).bundlename;
    idxTrl     = allSpks(su).idxTrlSing; % INDEXED TRIALS
    
    % DEFINE TIME WINDOW
    encTrig = allSpks(su).encTrigger(allSpks(su).hitsIdx, [1, 3]);
    encTrig = round(encTrig * 1000);
    
    % DIFFERENT TIME WINDOWS OF INTEREST
    switch  tw
        case 0                                                           % WHOLE TRIAL
            encTrig = [encTrig(:,1) encTrig(:,2)];
        case 1                                                           % PRE-CUE
            encTrig = [encTrig(:,1)-1000 encTrig(:,1)];
        case 2                                                           % PERI-CUE
            encTrig = [encTrig(:,1)-1000 encTrig(:,1)+1000];
        case 3                                                           % PERI-RESP
            encTrig = [encTrig(:,2)-1000 encTrig(:,2)+1000];
    end
    
    trlNum       = length(encTrig);
    
    numRppl      = zeros(1,trlNum);
    lenRppl      = cell(1,trlNum);
    rpplDen      = zeros(1,trlNum);
    
    
    % LOOP OVER TRIALS
    for trl = 1 : trlNum
        
        rppls = logical(allSpks(su).rppls{1,trl});
        rpplLen = allSpks(su).rpplLen{1,trl};
        if isempty(rppls)
            continue
        end
        
        %                     lenRpplAll = [];
        %                     trlRpplDen = zeros(1, [encTrig(trl,2)-encTrig(trl,1)+1]); % ripple density in that trial
        
        %                         if any(intersect([rppls(ripNum,1):rppls(ripNum,2)], [encTrig(trl,1) : encTrig(trl,2)]))                                     % is there any overlap with the current ripple and the current trial?
        
        numRppl(trl) = size(rpplLen,2) / size(rppls,2)*1000; % rppls per second                                                                           % one more ripple occurs within the current trial
        lenRppl{trl} = rpplLen;                              % length of the current ripple. is later averaged over all ripples in that trial
        rpplDen(trl) = mean(rppls);                          % ripple density
    end
    
    
    %% SPLIT TRIALS INTO IDX / NDX
    %  THEY ARE THEN SAVED IN VARIABLES LIKE "numRpplIdx". TRIAL IS
    %  A RANDOM EFFECT HERE. THIS PART IS SO THE VALUES FROM THE PREVIOUS STEP (WITHIN BUNDLE) ARE SAVED
    %  OVER BUNDLES
    
    %% INDEXED
    %  INDEXED TRIALS - NUMBER OF RIPPLES
    numRpplIdx = [numRpplIdx, numRppl( idxTrl )];
    
    %  INDEXED TRIALS - LENGTH OF RIPPLES
    temp = lenRppl( idxTrl );
    temp = temp(~cell2mat(cellfun(@isempty, temp, 'un', 0)));
    temp = [temp{:}];
    lenRpplIdx = [lenRpplIdx, temp];
    
    %  INDEXED TRIALS - RIPPLE DENSITY
    rpplDenIdx = [rpplDenIdx, rpplDen( idxTrl )];
    
    %             %  INDEXED TRIALS - SPIKE RIPPLE COINCIDENCE
    %             tmp = rpplSpkCoinc( idxTrl,:);
    %             coincIdx = [coincIdx; tmp(:)];
    
    %% NON-INDEXED
    %  NON-INDEXED TRIALS - NUMBER OF RIPPLES
    numRpplNdx = [numRpplNdx, numRppl(~idxTrl)];
    
    %  NON-INDEXED TRIALS - LENGTH OF RIPPLES 
    temp = lenRppl(~idxTrl);
    temp = temp(~cell2mat(cellfun(@isempty, temp, 'un', 0)));
    temp = [temp{:}];
    lenRpplNdx = [lenRpplNdx, temp];
    
    %  NON-INDEXED TRIALS - RIPPLE DENSITY
    rpplDenNdx = [rpplDenNdx, rpplDen(~idxTrl)];
    
    %             %  NON-INDEXED TRIALS - SPIKE RIPPLE COINCIDENCE
    %             tmp = rpplSpkCoinc(~idxTrl,:);
    %             coincNdx = [coincNdx; tmp(:)];
    
end % END OF SU LOOP

%% VISUALIZE
figure('units','normalized','outerposition',[0 0 1 1]);

% FIGURE TITLE
tits = '';
if tw == 0
    tits = strjoin({tits, 'Whole Trial |'});
elseif tw == 1
    tits = strjoin({tits, 'Pre-Cue |'});
elseif tw == 2
    tits = strjoin({tits, 'Peri-Cue |'});
elseif tw == 3
    tits = strjoin({tits, 'Peri-Resp |'});
end


tits = strjoin({tits, 'Only Microwire with maximum value'});
sgtitle(tits)

% RIPPLE NUMBER
subplot(211)
[p, ~] = ranksum(numRpplIdx, numRpplNdx, 'tail', 'right');
[p_perm] = perm_ranksum(numRpplIdx', numRpplNdx');
title({sprintf('Number of Ripples (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
    sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(numRpplIdx), mean(numRpplIdx), median(numRpplIdx)), ...
    sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(numRpplNdx), mean(numRpplNdx), median(numRpplNdx))});
hold on
hand1 = histogram(numRpplIdx, 0:1:10, 'Normalization','probability');
hand2 = histogram(numRpplNdx, 0:1:10, 'Normalization','probability');
hand2.FaceColor = [1,0,0];
plot([mean(numRpplIdx) mean(numRpplIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 4, 'linestyle', '-');
plot([mean(numRpplNdx) mean(numRpplNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '-');
plot([median(numRpplIdx) median(numRpplIdx)], get(gca, 'YLim'), 'color', 'b',     'linew', 4, 'linestyle', '--');
plot([median(numRpplNdx) median(numRpplNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '--');
legend('Indexed', 'Non-Indexed')
xlabel(tits);

% RIPPLE LENGTH
subplot(212)
[p, ~] = ranksum(lenRpplIdx, lenRpplNdx, 'tail', 'right');
[p_perm] = perm_ranksum(lenRpplIdx', lenRpplNdx');
title({sprintf('Length of Ripples (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
    sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(lenRpplIdx), nanmean(lenRpplIdx), nanmedian(lenRpplIdx)), ...
    sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(lenRpplNdx), nanmean(lenRpplNdx), nanmedian(lenRpplNdx))});

hold on
hand1 = histogram(lenRpplIdx, 30:5:120, 'Normalization','probability');
hand2 = histogram(lenRpplNdx, 30:5:120, 'Normalization','probability');
hand2.FaceColor = [1,0,0];
plot([nanmean(lenRpplIdx) nanmean(lenRpplIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 4, 'linestyle', '-');
plot([nanmean(lenRpplNdx) nanmean(lenRpplNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '-');
plot([nanmedian(lenRpplIdx) nanmedian(lenRpplIdx)], get(gca, 'YLim'), 'color', 'b',     'linew', 3, 'linestyle', '--');
plot([nanmedian(lenRpplNdx) nanmedian(lenRpplNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '--');
legend('Indexed', 'Non-Indexed')

% SAVE FIGURE
saveTits = '';
if tw == 0
    saveTits = strjoin({saveTits, 'WholeTrial_'});
elseif tw == 1
    saveTits = strjoin({saveTits, 'PreCue_'});
elseif tw == 2
    saveTits = strjoin({saveTits, 'PeriCue_'});
elseif tw == 3
    saveTits = strjoin({saveTits, 'PeriResp_'});
end

% if maxMW == 0
%     saveTits = strjoin({saveTits, 'allMW'});
% elseif maxMW == 1
%     saveTits = strjoin({saveTits, 'maxMW'});
% end

saveas(gcf,['\\analyse4.psy.gla.ac.uk\project0309\Luca\visu\rpplNum_rpplLen_', saveTits], 'png' )

% VISUALIZE RIPPLE DENSITY
close all
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(20,1,[4:20])

[p, ~] = ranksum(rpplDenIdx, rpplDenNdx, 'tail', 'right');
[p_perm] = perm_ranksum(lenRpplIdx', lenRpplNdx');
title({sprintf('Density of Ripples (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
    sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(rpplDenIdx), mean(rpplDenIdx), median(rpplDenIdx)), ...
    sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(rpplDenNdx), mean(rpplDenNdx), median(rpplDenNdx))});
hold on

hand1 = histogram(rpplDenIdx, 0:0.0025:0.1, 'Normalization','probability');
hand2 = histogram(rpplDenNdx, 0:0.0025:0.1, 'Normalization','probability');
hand2.FaceColor = [1,0,0];

plot([mean(rpplDenIdx) mean(rpplDenIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 2, 'linestyle', '-');
plot([mean(rpplDenIdx) mean(rpplDenIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 2, 'linestyle', '-');
plot([mean(rpplDenNdx) mean(rpplDenNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 2, 'linestyle', '-');

plot([median(rpplDenIdx) median(rpplDenIdx)], get(gca, 'YLim'), 'color', 'b',     'linew', 2, 'linestyle', '--');
plot([median(rpplDenNdx) median(rpplDenNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 2, 'linestyle', '--');

legend('Indexed', 'Non-Indexed')
%         sgtitle(tits)
xlabel(tits)

saveas(gcf,['\\analyse4.psy.gla.ac.uk\project0309\Luca\visu\rpplDen_', saveTits], 'png' )
% close all

%         % VISUALIZE RIPPLE COINCIDENCES
%         close all
%         figure('units','normalized','outerposition',[0 0 1 1]);
%
%         [p, ~] = ranksum(coincIdx, coincNdx, 'tail', 'right');
%         [p_perm] = perm_ranksum(coincIdx, coincNdx);
%         title({sprintf('Spike-Ripple Coincidences (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
%             sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(coincIdx), mean(coincIdx), median(coincIdx)), ...
%             sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(coincNdx), mean(coincNdx), median(coincNdx))});
%         hold on
%         hand1 = histogram(coincIdx, 0:1:10, 'Normalization','probability');
%         hand2 = histogram(coincNdx, 0:1:10, 'Normalization','probability');
%         hand2.FaceColor = [1,0,0];
%         plot([mean(coincIdx) mean(coincIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 4, 'linestyle', '-');
%         plot([mean(coincNdx) mean(coincNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '-');
%         plot([median(coincIdx) median(coincIdx)], get(gca, 'YLim'), 'color', 'b',     'linew', 4, 'linestyle', '--');
%         plot([median(coincNdx) median(coincNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '--');
%         legend('Indexed', 'Non-Indexed')
%         %         sgtitle(tits)
%         xlabel(tits);
% saveas(gcf,['\\analyse4.psy.gla.ac.uk\project0309\Luca\visu\rpplCoin_', saveTits], 'png' )

% close all