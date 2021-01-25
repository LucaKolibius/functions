%% RIPPLE STATISTICS WITH (1) RANDOM EFFECT AND AS A (2) WITHIN DIFFERENCE
%% (1) RANDOM EFFECT
whereAmI(0);
global prePath;

% load([prePath, 'Luca\data\allSbj\allSpksHZ_rppls.mat'], 'allSpks')
load([prePath, 'Luca\data\allSbj\allSpksHZ_rppls80to140.mat'], 'allSpks')
allBids = {allSpks.bidsID};
allBund = {allSpks.bundlename};
allSesh = {allSpks.sesh};
skpSpk  = zeros(1,length(allSpks));

%% PREALLOCATION
ndxNum = [];
ndxLen = [];
ndxDen = [];
idxNum = [];
idxLen = [];
idxDen = [];

for spk = 1 : length(allSpks)
    if sum(allSpks(spk).idxTrl) == 0
        continue
    end
    
    %% SKIP SPIKE
    if skpSpk(spk) == 1
        continue
    end
    
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
    bund    = allSpks(spk).bundlename;
    sameIdx = strcmp(bidsID, allBids) & strcmp(sesh, allSesh) & strcmp(bund, allBund);
    skpSpk(sameIdx) = 1;
    idxTrl  = any(vertcat(allSpks(sameIdx).idxTrl),1); % any trial that is indexed on that bundle
    
    %     encTrig = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1])*1000);
    
    for trl = 1:length(idxTrl)
        trlDur = allSpks(spk).encTrigger(allSpks(spk).hitsIdx(trl),3)-allSpks(spk).encTrigger(allSpks(spk).hitsIdx(trl),1);
        
        switch idxTrl(trl)
            case 0 % ndx trl
                ndxNum = [ndxNum, allSpks(spk).rpplNum(trl)/trlDur];
                ndxLen = [ndxLen, allSpks(spk).rpplLen{trl}];
                ndxDen = [ndxDen, mean(allSpks(spk).rpplDen{trl})];
            case 1 % idx trl
                idxNum = [idxNum, allSpks(spk).rpplNum(trl)/trlDur];
                idxLen = [idxLen, allSpks(spk).rpplLen{trl}];
                idxDen = [idxDen, mean(allSpks(spk).rpplDen{trl})];
        end
        
    end % OF TRL LOOP
    
end % OF SU LOOP

%%
%% (2) WITHIN DIFFERENCE
whereAmI(0);
global prePath;

load([prePath, 'Luca\data\allSbj\allSpksHZ_rppls.mat'], 'allSpks')
allBids = {allSpks.bidsID};
allBund = {allSpks.bundlename};
allSesh = {allSpks.sesh};

rpplNum.idx = [];
rpplNum.ndx = [];
rpplLen.idx = [];
rpplLen.ndx = [];
rpplDen.idx = [];
rpplDen.ndx = [];

skpSpk = zeros(1,length(allSpks));
for spk = 1 : length(allSpks)
    trlDur = [allSpks(spk).encTrigger(allSpks(spk).hitsIdx,3) - allSpks(spk).encTrigger(allSpks(spk).hitsIdx,1)]';
    
    allSpks(spk).rpplNum =     allSpks(spk).rpplNum ./trlDur;

    %% SKIP SPIKES ON THE SAME BUNLDE
    if skpSpk(spk) == 1
        continue
    end
    
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
    bund    = allSpks(spk).bundlename;
    sameIdx = strcmp(bidsID, allBids) & strcmp(sesh, allSesh) & strcmp(bund, allBund);
    skpSpk(sameIdx) = 1;
    idxTrl  = any(vertcat(allSpks(sameIdx).idxTrl),1); % any trial that is indexed on that bundle
    
    %% IF ANY TRIALS IN THAT BUNDLE ARE INDEXED
    if any(idxTrl)
        
        %% RIPPLE NUMBER
        rpplNum.idx  = [ rpplNum.idx  {allSpks(spk).rpplNum( idxTrl)}  ];
        rpplNum.ndx  = [ rpplNum.ndx  {allSpks(spk).rpplNum(~idxTrl)}  ];

        %% RIPPLE LENGTH
        rpplLen.idx  = [ rpplLen.idx  {allSpks(spk).rpplLen( idxTrl)}  ];
        rpplLen.ndx  = [ rpplLen.ndx  {allSpks(spk).rpplLen(~idxTrl)}  ];

        %% RIPPLE DENSITY
        rpplDen.idx  = [rpplDen.idx {cell2mat(cellfun(@mean, {allSpks(spk).rpplDen{ idxTrl}}, 'un', 0))}];
        rpplDen.ndx  = [rpplDen.ndx {cell2mat(cellfun(@mean, {allSpks(spk).rpplDen{~idxTrl}}, 'un', 0))}];
        
        
    end % IF THERE ARE INDEXED TRIALS
end % OF SU LOOP


% %% EMPIRICAL DIFFERENCE
% numComp = size(rpplDen.idx,2);
% diffLen_emp = [];
% diffNum_emp = [];
% diffDen_emp = [];
% 
% for comp = 1:numComp
%     %% LENGTH
%     empLenIdx = rpplLen.idx{comp};
%     empLenIdx = mean([empLenIdx{:}]);
%     
%     empLenNdx = rpplLen.ndx{comp};
%     empLenNdx = mean([empLenNdx{:}]);
%     
%     diffLen_emp = [diffLen_emp, empLenIdx - empLenNdx];
%     
%     %% NUMBER
%     empNumIdx = rpplNum.idx{comp};
%     empNumIdx = mean(empNumIdx);
%     
%     empNumNdx = rpplNum.ndx{comp};
%     empNumNdx = mean(empNumNdx);
%     
%     diffNum_emp = [diffNum_emp, empNumIdx - empNumNdx];
%     
%     %% DENSITY
%     empDenIdx = rpplDen.idx{comp};
%     empDenIdx = mean(empDenIdx);
%     
%     empDenNdx = rpplDen.ndx{comp};
%     empDenNdx = mean(empDenNdx);
%     
%     diffDen_emp = [diffDen_emp, empDenIdx - empDenNdx];
% end
% 
% %% PERMUTATION TEST FOR (2) WITHIN DIFFERENCE
% nperm = 10000;
% for perm = 1:nperm   
%     lineLength = fprintf('%d out of %d permutations done (%.2f%%).', perm, nperm, perm/nperm*100);
% 
%     lenDiffComp = [];
%     denDiffComp = [];
%     numDiffComp = [];
%     for comp = 1:numComp
%         %% LENGTH
%         idxLen = cell2mat(cellfun(@mean, rpplLen.idx{comp}, 'un', 0));
%         idxLen(isnan(idxLen)) = [];
%         ndxLen = cell2mat(cellfun(@mean, rpplLen.ndx{comp}, 'un', 0));
%         ndxLen(isnan(ndxLen)) = [];
%         
%         numIdx = size(idxLen,2);
%         
%         allLen = [idxLen ndxLen];
%         allLen = allLen(randperm(size(allLen,2)));
%         
%         idxLenPerm = mean(allLen(1:numIdx));
%         ndxLenPerm = mean(allLen(numIdx+1:end));
%         
%         lenDiffComp = [lenDiffComp idxLenPerm-ndxLenPerm];
% 
%         %% NUMBER
%         idxNum = rpplNum.idx{comp};
%         ndxNum = rpplNum.ndx{comp};
%         
%         numIdx = size(idxNum,2);
%         
%         allNum = [idxNum ndxNum];
%         allNum = allNum(randperm(size(allNum,2)));
% 
%         idxNumPerm = mean(allNum(1:numIdx));
%         ndxNumPerm = mean(allNum(numIdx+1:end));
%         
%         numDiffComp = [numDiffComp idxNumPerm-ndxNumPerm];
%         
%         %% DENSITY
%         idxDen = rpplDen.idx{comp};
%         ndxDen = rpplDen.ndx{comp};
%         
%         numIdx = size(idxDen,2);
%         
%         allDen = [idxDen ndxDen];
%         allDen = allDen(randperm(size(allDen,2)));
%         
%         idxDenPerm = mean(allDen(1:numIdx));
%         ndxDenPerm = mean(allDen(numIdx+1:end));
%         
%         denDiffComp = [denDiffComp, idxDenPerm-ndxDenPerm];
% 
%     end
%     numDiffPerm(perm) = mean(numDiffComp);
%     lenDiffPerm(perm) = nanmean(lenDiffComp);
%     denDiffPerm(perm) = mean(denDiffComp);
%             
%     
%     fprintf(repmat('\b',1,lineLength))
% 
% end
% %% NUMBER
% prctile(numDiffPerm, 95)
% mean(diffNum_emp)
% 
% sum(diffNum_emp > numDiffPerm)
% 
% %% LENGTH
% prctile(lenDiffPerm, 95)
% prctile(lenDiffPerm, 5)
% nanmean(diffLen_emp)
% 
% %% DENSITY
% prctile(denDiffPerm, 95)
% prctile(denDiffPerm,  5)
% mean(diffDen_emp)

%% VISUALISATION
figure('units','normalized','outerposition',[0 0 1 1]);

% RIPPLE NUMBER
subplot(211)
[p, ~] = ranksum(idxNum, ndxNum, 'tail', 'right');
[p_perm] = perm_ranksum(idxNum', ndxNum');
title({sprintf('Number of Ripples (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
    sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(idxNum), mean(idxNum), median(idxNum)), ...
    sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(ndxNum), mean(ndxNum), median(ndxNum))});
hold on
hand1 = histogram(idxNum, 0:1:10, 'Normalization','probability');
hand2 = histogram(ndxNum, 0:1:10, 'Normalization','probability');
hand2.FaceColor = [1,0,0];
plot([mean(idxNum) mean(idxNum)], get(gca, 'YLim'), 'color', 'b', 'linew', 4, 'linestyle', '-');
plot([mean(ndxNum) mean(ndxNum)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '-');
plot([median(idxNum) median(idxNum)], get(gca, 'YLim'), 'color', 'b',     'linew', 4, 'linestyle', '--');
plot([median(ndxNum) median(ndxNum)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '--');
legend('Indexed', 'Non-Indexed')

% RIPPLE LENGTH
subplot(212)
[p, ~] = ranksum(idxLen, ndxLen, 'tail', 'right');
[p_perm] = perm_ranksum(idxLen', ndxLen');
title({sprintf('Length of Ripples (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
    sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(idxLen), nanmean(idxLen), nanmedian(idxLen)), ...
    sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(ndxLen), nanmean(ndxLen), nanmedian(ndxLen))});

hold on
hand1 = histogram(idxLen, 30:5:120, 'Normalization','probability');
hand2 = histogram(ndxLen, 30:5:120, 'Normalization','probability');
hand2.FaceColor = [1,0,0];
plot([nanmean(idxLen) nanmean(idxLen)], get(gca, 'YLim'), 'color', 'b', 'linew', 4, 'linestyle', '-');
plot([nanmean(ndxLen) nanmean(ndxLen)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '-');
plot([nanmedian(idxLen) nanmedian(idxLen)], get(gca, 'YLim'), 'color', 'b',     'linew', 3, 'linestyle', '--');
plot([nanmedian(ndxLen) nanmedian(ndxLen)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '--');
legend('Indexed', 'Non-Indexed')


% VISUALIZE RIPPLE DENSITY
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(20,1,[4:20])

[p, ~] = ranksum(idxDen, ndxDen, 'tail', 'right');
[p_perm] = perm_ranksum(idxLen', ndxLen');
title({sprintf('Density of Ripples (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
    sprintf('Indexed Trials (#% d; M = %.2f; MED = %.2f)', length(idxDen), mean(idxDen), median(idxDen)), ...
    sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(ndxDen), mean(ndxDen), median(ndxDen))});
hold on

hand1 = histogram(idxDen, 0:0.0025:0.1, 'Normalization','probability');
hand2 = histogram(ndxDen, 0:0.0025:0.1, 'Normalization','probability');
hand2.FaceColor = [1,0,0];

% plot([mean(idxDen) mean(idxDen)], get(gca, 'YLim'), 'color', 'b', 'linew', 2, 'linestyle', '-');
% plot([mean(idxDen) mean(idxDen)], get(gca, 'YLim'), 'color', 'b', 'linew', 2, 'linestyle', '-');
% plot([mean(ndxDen) mean(ndxDen)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 2, 'linestyle', '-');
% 
% plot([median(idxDen) median(idxDen)], get(gca, 'YLim'), 'color', 'b',     'linew', 2, 'linestyle', '--');
% plot([median(ndxDen) median(ndxDen)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 2, 'linestyle', '--');

title(sprintf('Density of Ripples (p_{ranksum}) = %.2f', p))
yticks([0:0.05:0.2])
ylim([0 0.2])
legHand = legend('Indexed', 'Non-Indexed');
legHand.FontSize = 16;
legHand.FontWeight = 'bold';
xlabel('Density')
ylabel('Proportion');

figHand = gca;
figHand.FontSize = 16;
figHand.FontWeight = 'bold';