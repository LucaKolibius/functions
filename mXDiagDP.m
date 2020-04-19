function xDiagDP = mXDiagDP(encSpikeNumber_cueLocked_h, retSpikeNumber_cueLocked_h, hitMiss, su)

% create standardized SU rates in RSA orda
allTrials = size(encSpikeNumber_cueLocked_h,2);
if size(encSpikeNumber_cueLocked_h,1)>2
    [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_h, retSpikeNumber_cueLocked_h, hitMiss,0); % extracts the absolute spikenumbers for encoding and retrieval. CAREFUL: only hits and ordered by presentation time at encoding
else
    return
end

if size(tempCell_enc_h,2)>1 % in case there are no misses
    [~ , normENC, normRET] = normSpikeNumber(tempCell_enc_h, tempCell_ret_h); % normalizes spikenumber and correlates enc+ret for RSA
else
    return
end

%% number of FF/PP/FP for hits & misses
% hit
numFF_hit=0;
for i=1:size(encSpikeNumber_cueLocked_h,2)
    if strcmp(encSpikeNumber_cueLocked_h{1,i}, 'ff') && strcmp(encSpikeNumber_cueLocked_h{2,i}, 'hit')
        numFF_hit=numFF_hit+1;
    end
end

numPP_hit=0;
for i=1:size(encSpikeNumber_cueLocked_h,2)
    if strcmp(encSpikeNumber_cueLocked_h{1,i}, 'pp') && strcmp(encSpikeNumber_cueLocked_h{2,i}, 'hit')
        numPP_hit=numPP_hit+1;
    end
end
numPP_hit=numPP_hit+numFF_hit;

numFP_hit=0;
for i=1:size(encSpikeNumber_cueLocked_h,2)
    if strcmp(encSpikeNumber_cueLocked_h{1,i}, 'fp') && strcmp(encSpikeNumber_cueLocked_h{2,i}, 'hit')
        numFP_hit=numFP_hit+1;
    end
end
numFP_hit=numFP_hit+numPP_hit;
numHit=[numFF_hit numPP_hit numFP_hit];

% miss
numFF_miss=0;
for i=1:size(encSpikeNumber_cueLocked_h,2)
    if strcmp(encSpikeNumber_cueLocked_h{1,i}, 'ff') && strcmp(encSpikeNumber_cueLocked_h{2,i}, 'miss')
        numFF_miss=numFF_miss+1;
    end
end

numPP_miss=0;
for i=1:size(encSpikeNumber_cueLocked_h,2)
    if strcmp(encSpikeNumber_cueLocked_h{1,i}, 'pp') && strcmp(encSpikeNumber_cueLocked_h{2,i}, 'miss')
        numPP_miss=numPP_miss+1;
    end
end
numPP_miss = numPP_miss+numFF_miss;

numFP_miss = 0;
for i=1:size(encSpikeNumber_cueLocked_h,2)
    if strcmp(encSpikeNumber_cueLocked_h{1,i}, 'fp') && strcmp(encSpikeNumber_cueLocked_h{2,i}, 'miss')
        numFP_miss=numFP_miss+1;
    end
end
numFP_miss = numFP_miss+numPP_miss;
numMiss = [numFF_miss numPP_miss numFP_miss];

%%
encTS = normENC(su,:);
retTS = normRET(su,:);

% die WT dp wäre einfach nur encTS .* retTS
% hier will ich aber die WC raussuchen
% ff / pp / fp is the order
if strcmp(hitMiss,'hit')
    ffIdx = flip(1:numHit(1));
    ppIdx = flip(numHit(1)+1:numHit(2));
    fpIdx = flip(numHit(2)+1:numHit(3));
    
elseif strcmp(hitMiss,'miss')
    ffIdx = flip(1:numMiss(1));
    ppIdx = flip(numMiss(1)+1:numMiss(2));
    fpIdx = flip(numMiss(2)+1:numMiss(3));
end

retIdx = [ffIdx, ppIdx, fpIdx];
xDiagDP = encTS .* retTS(retIdx);

end % end of function