clear
close all
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_retResp.mat', 'allSpks');
nperm = 10000;
counter    = 0;
allIdxProx = [];
dt         = -5:1:5;
for su = 1:length(allSpks)
    idxTrl = [NaN(1,5) allSpks(su).idxTrl NaN(1,5)];
    if nansum(idxTrl)<2
        continue
    end
    counter = counter + 1;
    
    %% PERMUTATION TEST (WHY NOT IN THE BEGINNING FOR A TEST?)
    for perm = 1:nperm
        %         idxTrlPerm{counter}(perm,:) =
        tt = idxTrl(6:end-5);
        tt = tt(randperm(size(tt,2)));
        tt = [NaN(1,5) tt NaN(1,5)];
        numIdx  = find(tt==1);
        
        idxProx = [];
        for idx = 1:size(numIdx,2)
            idxProx = [idxProx; tt(numIdx(idx)-5:numIdx(idx)+5)];
        end
        idxProxPerm(counter, perm,:) = nanmean(idxProx,1);
    end
    
    %% EMPIRICAL PATTERN
    numIdx  = find(idxTrl==1);
    idxProx = [];
    for idx = 1:size(numIdx,2)
        idxProx = [idxProx; idxTrl(numIdx(idx)-5:numIdx(idx)+5)];
    end
    
%     figure(su); hold on;
%     for ii = 1:size(idxProx,1)
%         tt = idxProx(ii,:); tt(6) = NaN;
%         plot(dt, tt);
%     end
    idxProx    = nanmean(idxProx,1);
    allIdxProx = [allIdxProx; idxProx];
end
allIdxProx(:,6) = NaN;
reshape(idxProxPerm, [], 11)
idxProxPerm = reshape(idxProxPerm, [], 11);
idxProxPerm(:,6) = NaN;
proxTh = prctile(idxProxPerm,95,1);
allIdxProx(:,6) = NaN;
empProxProb = mean(allIdxProx,1);

%% pvals
% empProxProb
% idxProxPerm

mean(idxProxPerm > empProxProb, 1)

%% VISUALISATION
figure(1); clf;
plot(dt, empProxProb, 'linew', 2, 'color', [0 0 0])
hold on;

% ERROR IDX
shadeDT   = [dt, fliplr(dt)];
shade     = nanstd(allIdxProx,0,1) / sqrt(size(allIdxProx,1));
inBetween = [empProxProb-shade, fliplr(empProxProb+shade)];

inBetween([6 17]) = 0;
fillHand  = fill(shadeDT, inBetween, [0, 0, 0] );
fillHand.FaceAlpha = 0.5;
rectangle('Position', [-1 0.001 2 0.4], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1])


plot(dt, proxTh, '-', 'color', [1 0 0], 'linew', 3)
ylim([0 0.45])
yticks(0:0.1:0.5)
xlabel('Lag')
ylabel('Probability of reinstated trial')
title('Probability that trials next to an reinstated trial are also reinstated')
% mAx = gca;
% mAx.FontSize = 20;
% mAx.FontWeight = 'bold';
set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',25,'FontWeight','Bold', 'LineWidth', 2);

L(1) = plot(nan, nan, 'color', [0, 0, 0]);
L(2) = plot(nan, nan, 'color', [1, 0, 0]);
[legPos, hobj, ~, ~] = legend(L, {'Empirical proability', '95% threshold'}, 'FontSize',20, 'FontWeight', 'bold');
set(hobj,'LineWidth',15);
set(legPos, 'Position', [0.7608    0.8086    0.1962    0.1022])
legend('boxoff')
box off

