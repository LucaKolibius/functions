%% THIS IS CALLED IN lfp2rpplBund
function [idxTrlPow, ndxTrlPow, goodTrl, freqRes] = calcRipplPow (data, idxTrl, encTrig, favChan)

idxTrlPow = [];
ndxTrlPow = [];
goodTrl   = ones(1, length(idxTrl));

for trl = 1 : length(encTrig)
    
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
%     cfg.pad       = ceil(max(cellfun(@numel, data.time)/data.fsample));
    cfg.foi       = foi;
    cfg.trials    = trl;
    cfg.toi       = 'all';
    cfg.width     = 12;
    cfg.method    = 'wavelet';
    
    trlPow        = ft_freqanalysis(cfg, data);
    freqRes       = trlPow.freq;
    trlPow        = trlPow.powspctrm;        
    
    %% AR
    [isAR,~] = iqrAR(data.trial{trl},0);
    isAR     = squeeze(isAR(:,1:length(foi),:));
    trlPow(isAR) = nan;
    
    %% Chose Input Channel per Frequency
    temp = [];
    for freq = 1:6
        temp(freq,:) = trlPow(favChan(freq),freq,:);
    end
    temp                 = squeeze(temp);
    trlPow               = temp;
    
    %% Cut Wings
    trlPow(:,1:100)      = [];
    trlPow(:,end-99:end) = [];
    
    %% NORMALIZE TF PLOT
    baseline = trlPow(:,1:200);
    trlPow   = trlPow(:,201:end);
    
    meanBL = nanmean(baseline,2);
    stdBL  = nanstd(baseline,0,2);
    trlPow = (trlPow-meanBL)./stdBL;
    
    
    trlPow = nanmean(trlPow,2)';
    
    %     ONLY FOR POWERSPECTRUM
    %     adjPow = repmat(size(trlPowTime,2),1,8) - sum(isnan(trlPow),2); %
    %     trlPowTime = trlPowTime./adjPow; 
    
    switch idxTrl(trl)
        case 0
            ndxTrlPow = [ndxTrlPow; trlPow];
        case 1
            idxTrlPow = [idxTrlPow; trlPow];
    end
    
end % END OF TRIAL LOOP

end % END OF FUNCTION