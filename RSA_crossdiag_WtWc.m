clear all
close all
cd X:\Luca\data\analysis
load('SQ.mat');
fetch_allSubj
SPSSinput= zeros(size(allSubj,1),4); % WT-hit WC-hit WT-miss WC-miss
for sbj = 1:size(allSubj,1)
try
    cd X:/Luca/data
catch
    cd /media/ldk898/rds-share/Luca/data
end
subjID = allSubj{sbj};
mSubject = subjID(1:end-3);
mSession = subjID(end-1:end);

% if the session name is called 1b then this line prevents an error during cd
mSubject(regexp(mSubject,'_')) = [];
if isempty(regexp(mSession,'S', 'ONCE'))
    mSession = ['S', mSession];
end

cd(mSubject)
cd(mSession)
abc = dir;
cd(abc(3).name)
cd advancedAnalysis\RSA\oldSize\hippocampus
abc = dir('RSAmat_hipp*');
if isempty(abc) % this means we don't have a hippocampus RSA / wires are all in cortex
    continue
end
load(abc.name, 'RSAhits_CueCue_h', 'RSAmiss_CueCue_h');

if size(RSAmiss_CueCue_h,1)<4 % at least 3 miss trials to be considered
    continue
end

for sqNum = [1 5 9]
    % hit square 1,5,9
    sqHit = RSAhits_CueCue_h(SQ(sbj).mask_hits==sqNum); % extracts empirical values for square sq out of RSA based on mask
    sqHit = reshape(sqHit,SQ(sbj).size_hits(:,sqNum)'); % resize to original size
    
    % miss square 1,5,9
    sqMis = RSAmiss_CueCue_h(SQ(sbj).mask_miss==sqNum);
    sqMis = reshape(sqMis,SQ(sbj).size_miss(:,sqNum)');
    
    
    SPSSinput(sbj, 1)  =  SPSSinput(sbj, 1) +  1/3 * mean(sqHit(logical(eye(size(sqHit)))));
    SPSSinput(sbj, 2)  =  SPSSinput(sbj, 2) +  1/3 * mean(sqHit(logical(rot90(eye(size(sqHit))))));
    SPSSinput(sbj, 3)  =  SPSSinput(sbj, 1) +  1/3 * mean(sqMis(logical(eye(size(sqMis)))));
    SPSSinput(sbj, 4)  =  SPSSinput(sbj, 2) +  1/3 * mean(sqMis(logical(rot90(eye(size(sqMis))))));
end
end

idx = SPSSinput == 0;
idx = sum(idx,2) == 4;
SPSSinput(idx,:) = [];