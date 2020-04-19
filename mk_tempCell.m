% extracts all the number of spikes of hits or misses. First extracts 'ff'
% trials, then 'pp' and 'fp'

% example for cueLocked & hits:
% mk_tempCell(allTrials, encSpikeNumber_cueLocked, retSpikeNumber_cueLocked, hitOrMiss, myFlag)
function [tempCell_enc, tempCell_ret]=mk_tempCell(allTrials, encSpikeNumber, retSpikeNumber, hitOrMiss, myFlag)
% encoding
tempCell_enc={};
for i=1:allTrials
    if strcmp(encSpikeNumber{1,i},'ff') && strcmp(encSpikeNumber{2,i}, hitOrMiss)
        tempCell_enc(:,size(tempCell_enc,2)+1)=encSpikeNumber{3:end,i};
    end
end

for i=1:allTrials
    if strcmp(encSpikeNumber{1,i},'pp') && strcmp(encSpikeNumber{2,i}, hitOrMiss)
        tempCell_enc(:,size(tempCell_enc,2)+1)=encSpikeNumber{3:end,i};
    end
end

for i=1:allTrials
    if strcmp(encSpikeNumber{1,i},'fp') && strcmp(encSpikeNumber{2,i}, hitOrMiss)
        tempCell_enc(:,size(tempCell_enc,2)+1)=encSpikeNumber{3:end,i};
    end
end

% retrieval
tempCell_ret={};
for i=1:allTrials
    if strcmp(retSpikeNumber{1,i},'ff') && strcmp(retSpikeNumber{2,i}, hitOrMiss)
        tempCell_ret(:,size(tempCell_ret,2)+1)=retSpikeNumber{3:end,i};
    end
end

for i=1:allTrials
    if strcmp(retSpikeNumber{1,i},'pp') && strcmp(retSpikeNumber{2,i}, hitOrMiss)
        tempCell_ret(:,size(tempCell_ret,2)+1)=retSpikeNumber{3:end,i};
    end
end

for i=1:allTrials
    if strcmp(retSpikeNumber{1,i},'fp') && strcmp(retSpikeNumber{2,i}, hitOrMiss)
        tempCell_ret(:,size(tempCell_ret,2)+1)=retSpikeNumber{3:end,i};
    end
end

% delete flagged rows
if myFlag~=0
tempCell_enc(myFlag==1,:)=[];
tempCell_ret(myFlag==1,:)=[];
end
end
