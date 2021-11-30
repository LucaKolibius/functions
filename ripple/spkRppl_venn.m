trlLen = 20000;
coinc_encIdx = [];
coinc_encNdx = [];

for su = 1:length(allSpks)
    
    if allSpks(su).iu == 0
        continue
    end
    
    rips    = allSpks(su).rpplTms;
    spks    = allSpks(su).spks/1000;
    trig    = allSpks(su).encTrigger(allSpks(su).hitsIdx,:);
    spksSeg       = insertSpiketimes2(trig, spks, [1 1], [0 20])';
    spksSeg       = cellfun(@(x) x*1000,   spksSeg, 'un', 0);
    spksSeg       = cellfun(@(x) round(x), spksSeg, 'un', 0);
    idxT = allSpks(su).idxTrl;
    
    k = 1;
    while isempty(rips)
        rips = allSpks(su-k).rpplTms;
        k = k+1;
    end
    
    coinc_norm = [];
    for trl = 1:size(idxT,2)
        spkTrl = spksSeg{trl,1};
        spkTrl(spkTrl == 0) = 1; % 0 due to rounding
        
        % IMAGINE A VENN DIAGRAM
        p_spk = size(spksSeg{trl,1},1) / trlLen;
        p_rip = mean(rips(trl,:),2);
        coinc = sum(rips(trl,spkTrl)) / trlLen;
        
        coinc_norm(trl,1) = coinc / (p_spk + p_rip - coinc);
        
        if p_spk == 0 || p_rip == 0
            coinc_norm(trl,1) = NaN;
        end
        
    end
    
    coinc_encIdx = [coinc_encIdx; coinc_norm( idxT)];
    coinc_encNdx = [coinc_encNdx; coinc_norm(~idxT)];
end
figure(1); clf;
dt = 0:0.005:0.1;
subplot(211); histogram(coinc_encIdx, dt, 'Normalization', 'Probability');
subplot(212); histogram(coinc_encNdx, dt, 'Normalization', 'Probability');

ranksum(coinc_encIdx, coinc_encNdx)

nanmedian(coinc_encIdx)
nanmedian(coinc_encNdx)

%% ripples missing for su601 & 605