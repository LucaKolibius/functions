% THIS ONE PRESERVES THE ENCODING SEQUENCE (RETRIEVAL SEQUENCE IS
% AUTOMATICALLY MATCHED TO ENCODING SEQUENCE THROUGH "LOADLOGS"
% extracts all the number of spikes of hits or misses.

% example for cueLocked & hits:
% mk_tempCell2(allTrials, encSpikeNumber_cueLocked, retSpikeNumber_cueLocked, hitOrMiss)
function [tempCell_enc, tempCell_ret]=mk_tempCell2(allTrials, encSpikeNumber, retSpikeNumber, hitOrMiss)
% encoding
tempCell_enc={};
for i=1:allTrials
    if strcmp(encSpikeNumber{2,i}, hitOrMiss)
        tempCell_enc(:,size(tempCell_enc,2)+1)=encSpikeNumber{3:end,i};
    end
end

% retrieval
tempCell_ret={};
for i=1:allTrials
    if strcmp(retSpikeNumber{2,i}, hitOrMiss)
        tempCell_ret(:,size(tempCell_ret,2)+1)=retSpikeNumber{3:end,i};
    end
end
end
