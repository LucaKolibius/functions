% %% just looking at trial order
% clear
% load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_retResp.mat', 'allSpks');
% 
% nperm = 10000;
% counter = 0;
% pPermHi = [];
% pPermLo = [];
% for su = 1:length(allSpks)
%     
% su
% 
% idxTrl = allSpks(su).idxTrl;
% if sum(idxTrl) < 2
%     continue
% end
% counter = counter + 1;
% diffIdx = find(idxTrl);
% diffIdx = diffIdx(2:end) - diffIdx(1:end-1);
% diffIdx = mean(diffIdx);
% 
% diffIdxPerm = zeros(1,nperm);
% for perm = 1:nperm
%     idxTrlPerm = idxTrl(randperm(size(idxTrl,2)));
%     
%     tt                = find(idxTrlPerm);
%     tt                = tt(2:end) - tt(1:end-1);
%     diffIdxPerm(perm) = mean(tt);
%     
% end
% 
% pPerm(counter) = 1 - (sum(diffIdx<diffIdxPerm) / nperm);
% 
% end
% 
% figure(1); clf;
% plot(pPerm, 'o-', 'linew', 2); hold on;
% xlim([1 length(pPerm)])
% plot(get(gca, 'Xlim'), [0.05 0.05], '--', 'color', [1 0 0], 'linew', 2)
% % plot(get(gca, 'Xlim'), [0.975 0.975], '--', 'color', [1 0 0], 'linew', 2)
% xlabel('IU')
% title({sprintf('IU with indexed trials that are significantly closer together (#%d of #%d)',sum(pPerm < 0.05), length(pPerm)), 'Order of Trials'})
% ylabel('p-Value (H1: empirical distance is lower than permuted)')

%% looking at actual time
clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_retResp.mat', 'allSpks');

onlyHits = 1;
nperm    = 10000;
counter  = 0;
pPermHi  = [];
pPermLo  = [];
for su = 1:length(allSpks)
    
su

idxTrl = allSpks(su).idxTrl;
if sum(idxTrl) < 2
    continue
end
counter = counter + 1;

switch onlyHits
    case 0
        encTrl = allSpks(su).encTrigger(:,1);
    case 1
        encTrl = allSpks(su).encTrigger(allSpks(su).hitsIdx,1);
end

idx     = find(idxTrl);
tt = encTrl(idx);
tt = tt(2:end) - tt(1:end-1);
idxTime(counter) = mean(tt);

idxTimePerm{counter} = zeros(1,nperm);
for perm = 1:nperm
    idxTrlPerm = idxTrl(randperm(size(idxTrl,2)));
    
    idx                        = find(idxTrlPerm);
    tt                         = encTrl(idx);
    tt                         = tt(2:end) - tt(1:end-1);
    idxTimePerm{counter}(perm) = mean(tt);
 
end

pPerm(counter) = 1 - (sum(idxTime(counter)<idxTimePerm{counter}) / nperm);
th(counter)    = prctile(idxTimePerm{counter},5); % needs to be under that
end

sum(idxTime<=th);

%% CHECK FOR SIGNIFICANCE
for perm = 1:nperm
    expectedNull(perm) = 0; 
    for iu = 1:counter
        idx = randi(nperm); % grab random permutation
        
        if idxTimePerm{iu}(idx) <= th(iu)
            expectedNull(perm) = expectedNull(perm) + 1;
        end
        
    end
end

prctile(expectedNull,95)

%% VISU
figure(2); clf;
plot(pPerm, 'o-', 'linew', 2, 'color', [0 0 0]); hold on;
xlim([1 length(pPerm)])
plot(get(gca, 'Xlim'), [0.05 0.05], '-', 'color', [1 0 0], 'linew', 2)
% plot(get(gca, 'Xlim'), [0.975 0.975], '--', 'color', [1 0 0], 'linew', 2)
xlabel('Index Unit')
title({sprintf('IU with indexed trials that are significantly closer together (#%d of #%d)',sum(pPerm < 0.05), length(pPerm)), 'in Time [s] (binomial test p = 0.085)'})
ylabel({'p-Value (H1: empirical distance', 'is significantly lower than chance)'})

xlim([1 21])
mAx = gca;
mAx.FontSize = 20;
mAx.FontWeight = 'bold';
