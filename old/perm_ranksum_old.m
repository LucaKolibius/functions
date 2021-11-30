% behaves weirdly
function p = perm_ranksum_old(vec1, vec2)
nperm   = 10000;
numVec1 = size(vec1,1);
numVec2 = size(vec2,1);

%% EMPIRICAL RANKSUM
vec1 = [vec1, ones(size(vec1,1),1)*1]; 
vec2 = [vec2, ones(size(vec2,1),1)*2]; 

empDat = [vec1; vec2];

empDat(:,1) = floor(tiedrank([empDat(:,1)]));

rank1 = sum(empDat(empDat(:,2) == 1,1)) / numVec1;
rank2 = sum(empDat(empDat(:,2) == 2,1)) / numVec2;

relRank = rank1/rank2;

%% PERMUTATION TEST
relRankperm = zeros(1,nperm);
for perm = 1 : nperm
permDx  = randperm(numVec1 + numVec2);
permDat = [empDat(permDx,1), empDat(:,2)];

permDat(:,1) = floor(tiedrank([permDat(:,1)]));

permRank1 = sum(permDat(permDat(:,2) == 1,1)) / numVec1;
permRank2 = sum(permDat(permDat(:,2) == 2,1)) / numVec2;

relRankperm(perm) = permRank1 / permRank2;

end

%% ESTIMATE P-VALUE
p = 1 - (sum(relRank >= relRankperm) / nperm);

end  % END OF FUNCTION
