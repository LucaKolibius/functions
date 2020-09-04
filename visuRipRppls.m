load('X:\Luca\data\allSbj\rpplBund.mat')
load('X:\Luca\data\allSbj\allSpks.mat', 'allSpks');
lfpFold = 'X:\Luca\data\microLFP\';
tw = 500;
for it = 1 : length(rpplBund)
    
    %% LOAD RIPPLE DATA
    bidsID   = rpplBund(it).bidsID;
    sesh     = rpplBund(it).sesh;
    bundname = rpplBund(it).bundname;
    rppls    = rpplBund(it).rppls;
    
    %% LOAD SPIKETIMES FROM THAT SESSION
    idx = and(strcmp({allSpks.bidsID}, bidsID), strcmp({allSpks.sesh}, sesh) ); % index
    spks = vertcat(allSpks(idx).spks); % concatenate
    spks = sort(spks); % sort in time
    
    %% LOAD LFP DATA
    load([lfpFold, bidsID, '_', sesh, '_onlyMicroLFP_RAW_1000DS_SPKINT.mat'], 'data')
    
    % PICK THE RELEVANT MW FROM THE LFP
    allMW = cellfun(@(x) x(1:end-1), data.label, 'un', 0);
    idx   = strcmp(allMW, bundname);
    
    cfg             = [];
    cfg.channel     = data.label(idx);         % select MW from current bundle
    curMW           = ft_selectdata(cfg, data);
    
    cfg = [];
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [80 140];
    curMW_bp        = ft_preprocessing(cfg,curMW);
    
    mw = 2;  % MW on which the ripples were found (1-8)
    rpplMW = rppls{mw};

    figure('units','normalized','outerposition',[0 0 1 1]);

    for rip = 1: length(rpplMW)
        clf
    subplot(211)
    plot(curMW.trial{1}(mw, rpplMW(rip,1)-tw : rpplMW(rip,2)+tw), 'linew', 2); hold on
    plot([tw tw], get(gca, 'Ylim'), 'color', 'k', 'linew', 3);
    plot([rpplMW(rip,2)-rpplMW(rip,1)+tw rpplMW(rip,2)-rpplMW(rip,1)+tw;], get(gca, 'Ylim'), 'color', 'k', 'linew', 3);
    
    %% ADD SOME SPICY SPKS
    spksTW = spks(spks >= rpplMW(rip,1)-tw & spks <= rpplMW(rip,2)+tw); % spikes in that time window
    spksTW = spksTW - (rpplMW(rip,1)-tw) + 1;
    
    for ss = 1 : length(spksTW)
            plot([spksTW spksTW], get(gca, 'YLim'), 'color', 'r', 'linewidth', 0.1);
    end
    hold off

    %%
    subplot(212)
    plot(curMW_bp.trial{1}(mw, rpplMW(rip,1)-tw : rpplMW(rip,2)+tw), 'linew', 2); hold on
    plot([tw tw], get(gca, 'Ylim'), 'color', 'k', 'linew', 3);
    plot([rpplMW(rip,2)-rpplMW(rip,1)+tw rpplMW(rip,2)-rpplMW(rip,1)+tw;], get(gca, 'Ylim'), 'color', 'k', 'linew', 3);
    
    %% ADD SPIKES
     for ss = 1 : length(spksTW)
            plot([spksTW spksTW], get(gca, 'YLim'), 'color', 'r', 'linewidth', 0.1);
     end
    
    hold off
 
    saveas(gca, ['X:\Luca\rppl_spkint\rppl_spkint', num2str(rip)], 'png');
%     pause(2)
    end
    
end


end

