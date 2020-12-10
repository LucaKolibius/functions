%% THIS IS CALLED IN lfp2rpplBund
function [idxTrlPow, ndxTrlPow, goodTrl, freqRes] = calcRipplPow (data, idxTrl, encTrig, favChan)

idxTrlPow = [];
ndxTrlPow = [];
goodTrl   = ones(1, length(idxTrl));

for trl = 1 : length(encTrig)
    
    %% AR
    [isAR,~] = iqrAR(data.trial{trl},0);
    isAR = squeeze(isAR(1,1,:));
    
    %     if any(isAR)
    %         goodTrl(trl) = 0;
    %         continue
    %     end
    
    % FIELDTRIP POWER
    %     cfg.method    = 'mtmfft';
    %     cfg.taper     = 'hanning';
    
    foi = logspace(log10(1),log10(140));
    foi = foi(45:50);
    cfg = [];
    cfg.output    = 'pow';
    cfg.channel   = 'all';
    cfg.pad       = ceil(max(cellfun(@numel, data.time)/data.fsample));
    cfg.foi       = foi;
    cfg.trials    = trl;
    cfg.toi       = 'all';
    cfg.width     = 12;
    cfg.method    = 'wavelet';
    
    trlPow        = ft_freqanalysis(cfg, data);
    freqRes       = trlPow.freq;
    trlPow        = trlPow.powspctrm;
    
    
    
    trlPow        = trlPow/size(data.time,2); % adapt power spectrum by number of samples
    
    

    
    temp = [];
    for freq = 1:6
        temp(freq) = trlPow(favChan(freq),freq);
    end
    trlPow = temp;
    
    switch idxTrl(trl)
        case 0
            ndxTrlPow = [ndxTrlPow; trlPow];
        case 1
            idxTrlPow = [idxTrlPow; trlPow];
    end
    
end % END OF TRIAL LOOP

end % END OF FUNCTION