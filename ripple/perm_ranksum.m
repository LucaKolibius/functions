% behaves weirdly
function p = perm_ranksum(vec1, vec2)

%% MEDIAN
% nperm   = 10000;
% numVec1 = size(vec1,1);
% numVec2 = size(vec2,1);
% 
% empDiff = median(vec1) - median(vec2);
% 
% allDat = [vec1; vec2];
% for perm = 1 : nperm
%     randIdx = randperm(length(allDat));
%     allDat = allDat(randIdx);
%     
%     v1perm = allDat(1:numVec1);
%     v2perm = allDat(numVec1+1:end);
%     
%     permDiff(perm) = median(vec1) - median(vec2);
% end
% 
% % p = 1 - mean(empDiff >= permDiff) ; % I think that would work just fine
% p = 1 - (sum(empDiff >= permDiff) / nperm);


%% MEAN
nperm   = 10000;
numVec1 = size(vec1,1);
numVec2 = size(vec2,1);

empDiff = mean(vec1) - mean(vec2);

allDat = [vec1; vec2];
for perm = 1 : nperm
    randIdx = randperm(length(allDat));
    allDat = allDat(randIdx);
    
    v1perm = allDat(1:numVec1);
    v2perm = allDat(numVec1+1:end);
    
    permDiff(perm) = mean(vec1) - mean(vec2);
end

% p = 1 - mean(empDiff >= permDiff) ; % I think that would work just fine
p = 1 - (sum(empDiff >= permDiff) / nperm);

end