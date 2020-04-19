datapath = 'X:\Luca\data\allSbj';
dat = dir([datapath, filesep, 'allSU_*']);

numSU = [];
for ix = 1:length(dat)
load([dat(ix).folder, filesep, dat(ix).name])
countSU = size(allSU,1);
mName = dat(ix).name;
mName([1:6 end-3:end]) = [];
numSU = [numSU, [{mName}; countSU]];
end

sum([numSU{2,:}]) % 403 SU
sum([numSU{3,:}]) % 13 IU

backup = numSU;
% numSU(4,:) = num2cell(cellfun(@(x,y) y/x, numSU(2,:), numSU(3,:), 'UniformOutput', false))
numSU(3,cellfun(@isempty, numSU(3,:))) = {0};

numSU(4,:) = num2cell([numSU{3,:}] ./ [numSU{2,:}]);

nanmean([numSU{4,cell2mat(numSU(4,:))~= 0}])
sum([numSU{2,cell2mat(numSU(4,:))~= 0}])

%%

allWind = dir('X:\Luca\indexSUvisu\varWind*');
allStats = [];

for fold = 1:size(allWind,1)
    load([allWind(fold).folder, filesep, allWind(fold).name, filesep, 'globalTest'])
%     if sumnum>9500
        encWind = allWind(fold).name; retWind = allWind(fold).name;
        encWind = encWind( regexp(encWind, 'ENC_')+4: regexp(encWind, '_RET_')-1 );
        retWind = retWind( regexp(retWind, 'RET_')+4: regexp(retWind, '_DP_' )-1 );
        
        newStats = [sumnum; empNumFP2; th_fp2_glob; {encWind}; {retWind}];
        allStats = [allStats, newStats];
%     end
end

%% labels
labels = {};
for i = 1:27
[timeWindowENC, ~, lockingENC, ~, ~] = findIndexSU_timeWindows(i, 1);
labels = [labels {num2str([timeWindowENC, lockingENC])}];
end

figure;
subplot(121)
allSumNum = [allStats{1,:}];
% allSumNum(allSumNum>10000) = allSumNum(allSumNum>10000)/10;
% allSumNum = 1-(allSumNum/10000);
% [pthr, pcor, padj] = fdr(allSumNum)
allSumNum = reshape(allSumNum, 27, 27)';
imagesc(allSumNum);
axis square
xticks(1:1:27)
xticklabels(labels);
xtickangle(90)
yticks(1:1:27)
yticklabels(labels)

subplot(122)
% [mRow, mColumn] = find(allSumNum < 9500);
% allSumNum(mRow, mColumn) = 0;
allSumNum = [allStats{1,:}];
allSumNum(allSumNum<9500) = 0;
allSumNum = reshape(allSumNum, 27, 27)';
imagesc(allSumNum);
axis square
xlabel('Retrieval Window')
ylabel('Encoding Window')
xticks(1:1:27)
xticklabels(labels);
xtickangle(90)
yticks(1:1:27)
yticklabels(labels)

%%
for fold = 1:size(allWind,1)
    load([allWind(fold).folder, filesep, allWind(fold).name, filesep, 'globalTest'])
    if sumnum>=9900
        visuSU([allWind(fold).folder, filesep, allWind(fold).name])
    end
end

%%
allWind = dir('X:\Luca\indexSUvisu\varWind*');

perms_th = [];
for fold = 1:size(allWind,1)
    load([allWind(fold).folder, filesep, allWind(fold).name, filesep, 'globalTest'])
    thresh = empNumFP2*ones(size(permNumFP2_glob));
    perms_th = [perms_th, [permNumFP2_glob; thresh]];
end

nperm = 10000;
signTWperm_count = zeros(1,nperm);
numTW = size(allWind,1);
signTWperm = perms_th(1,:) <= perms_th(2,:);
for perm = 1:nperm
    tic
    randIdx = randperm(size(perms_th,2)); % random permutation index
    signTWperm = signTWperm(randIdx); % shuffle
    temp = signTWperm(1:10000*numTW); % take the first nperm*numTW number of shuffled entries
    temp = reshape(temp, 10000,numTW); % resize to make summing quicker
    issign = sum(temp, 1); % number of significant permutated cluster
    signTWperm_count(perm) = sum(issign>=9500); % count how many TW under the null become significant
    toc
end

