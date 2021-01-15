%% THIS IS CALLED IN lfp2rpplBund
function [idxTrlPow, ndxTrlPow, freqRes] = calcRipplPow (data, dataEx, idxTrl, encTrig, favChan)

idxTrlPow = [];
ndxTrlPow = [];

for trl = 1 : length(encTrig)
    
    %     if any(isAR)
    %         goodTrl(trl) = 0;
    %         continue
    %     end
    
    % FIELDTRIP POWER
    %     cfg.method    = 'mtmfft';
    %     cfg.taper     = 'hanning';
    
    cfg = [];
    cfg.output    = 'pow';
    cfg.channel   = 'all';
    cfg.pad       = ceil(max(cellfun(@numel, dataEx.time)/dataEx.fsample));
    cfg.foi       = 70:1:150;
    cfg.trials    = trl;
    cfg.toi       = 'all';
    cfg.width     = 12;
    cfg.method    = 'wavelet';
    
    trlPow        = ft_freqanalysis(cfg, dataEx);
    freqRes       = trlPow.freq;
    trlPow        = trlPow.powspctrm;        
    
    %% SNIP WINGS
    trlPow(:, :, 1:100)      = [];
    trlPow(:, :, end-99:end) = [];
    
    %% AR
    [isAR,~]     = iqrAR(data.trial{trl},0);
    isAR         = squeeze(isAR(:,1:length(cfg.foi),:));
    trlPow(isAR) = nan;
    
%     %% Chose Input Channel per Frequency
%     temp = [];
%     for freq = 1:6
%         temp(freq,:) = trlPow(favChan(freq),freq,:);
%     end
%     temp                 = squeeze(temp);
%     trlPow               = temp;
%     
    %% CHOSE INPUT CHANNEL
    trlPow = squeeze(trlPow(favChan,:,:));
    trlPow = nanmean(trlPow,2);
%     %% NORMALIZE TF PLOT (OLD)
%     baseline = trlPow(:,1:200);
%     trlPow   = trlPow(:,201:end);
%     
%     meanBL = nanmean(baseline,2);
%     stdBL  = nanstd(baseline,0,2);
%     trlPow = (trlPow-meanBL)./stdBL;
%     
%     
%     trlPow = nanmean(trlPow,2)';

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