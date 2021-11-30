function orda = reorderTrials(encSpikeNumber_tw_loc)
% this function reorders the trial according to their presentation sequence
% the first column is the presentation position in the whole experiment
% the second column shows the position after reordering into hits (ff-> pp -> fp) => miss (ff-> pp -> fp)
%   this is the order we get from "mk_tempCell"
% the third column is the category (ff/pp/fp)
% the fourth column is hits (1) or miss (0)

missIdx = ~cellfun(@isempty, strfind(encSpikeNumber_tw_loc{2,:}, 'miss'));
hitsIdx = ~cellfun(@isempty, strfind(encSpikeNumber_tw_loc{2,:}, 'hit'));

ffIdx = ~cellfun(@isempty, strfind(encSpikeNumber_tw_loc{1,:}, 'ff'));
ppIdx = ~cellfun(@isempty, strfind(encSpikeNumber_tw_loc{1,:}, 'pp'));
fpIdx = ~cellfun(@isempty, strfind(encSpikeNumber_tw_loc{1,:}, 'fp'));
trialNum = size(encSpikeNumber_tw_loc,2);

orda = [1:trialNum]'; 
posi = 0;
for i = 1:trialNum
    if hitsIdx(i) == 1 & ffIdx(i) == 1
        posi = posi +1;
        orda(i,2) = posi;
                orda(i,3) = 1;
        orda(i,4) = 1;
    end
end

for i = 1:trialNum
    if hitsIdx(i) == 1 & ppIdx(i) == 1
         posi = posi +1;
        orda(i,2) = posi;
                orda(i,3) = 2;
        orda(i,4) = 1;
    end
end

for i = 1:trialNum
    if hitsIdx(i) == 1 & fpIdx(i) == 1
         posi = posi +1;
        orda(i,2) = posi;
                orda(i,3) = 3;
        orda(i,4) = 1;
    end
end

for i = 1:trialNum
    if missIdx(i) == 1 & ffIdx(i) == 1
         posi = posi +1;
        orda(i,2) = posi;
        orda(i,3) = 1;
        orda(i,4) = 0;
    end
end

for i = 1:trialNum
    if missIdx(i) == 1 & ppIdx(i) == 1
         posi = posi +1;
        orda(i,2) = posi;
                orda(i,3) = 2;
        orda(i,4) = 0;
    end
end

for i = 1:trialNum
    if missIdx(i) == 1 & fpIdx(i) == 1
         posi = posi +1;
        orda(i,2) = posi;
                orda(i,3) = 3;
        orda(i,4) = 0;
    end
end

orda = sortrows(orda,2);