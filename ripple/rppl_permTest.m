% %% test quantity
% idx = [IUidx; GUidx];
% ndx = [IUndx; GUndx; SUnbd];
% [p, h] = ranksum(idx,ndx, 'tail', 'right') 
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% hold on
% histogram(idx, 0:3:100, 'Normalization','probability')
% histogram(ndx, 0:3:100, 'Normalization','probability')
% % histogram(idx, 0:1:30, 'Normalization','probability')
% % histogram(ndx, 0:1:30, 'Normalization','probability')
% plot([nanmedian(idx), nanmedian(idx)], get(gca,'YLim'), 'linew', 3, 'color', 'b');
% plot([nanmedian(ndx), nanmedian(ndx)], get(gca,'YLim'), 'linew', 3, 'color', 'r');
% 
% nidx   = size(idx,1);
% allDat = [idx; ndx]; 
% nperm  = 10000;
% medDif = zeros(1,nperm);
% for perm = 1:nperm
%     
%     randIdx = randperm(size(allDat,1));
%     allDat = allDat(randIdx);
%     
%     idxPerm = allDat(1:size(idx,1));
%     ndxPerm = allDat(size(idx,1)+1:end);
%     
%     medDif(perm) = nanmedian(idxPerm)-nanmedian(ndxPerm);
%     
% end
% 
% thUp = prctile(medDif,95);
% empDiff = nanmedian(idx)-nanmedian(ndx);
% 1-sum(empDiff>medDif)/nperm
% figure('units','normalized','outerposition',[0 0 1 1])
% histogram(medDif)
% hold on
% plot([empDiff empDiff], get(gca, 'YLim'), 'linew', 3, 'color', 'r');