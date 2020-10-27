%% SPIKE PPC 
%  This does not take any artefacts into account (they have to be NaN'ed before
clear
server = 0;

switch server
    case 0     % RUN ON MY OWN LAPTOP
        addpath('X:\Common\toolboxes\fieldtrip-20200310');
        addpath('X:\Luca\toolboxes\CircStat2012a');
        ft_defaults
        %         allMicro = 'X:\Luca\data\microLFP\'; % non interpolated LFP
        allMicro = 'X:\Luca\data\microLFP_dmeanOrth\'; % demeaned and orth. LFP
        load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks');
        
            case 1         % RUN ON THE SERVER

                addpath('/castles/nr/projects/h/hanslmas-ieeg-compute/Common/toolboxes/fieldtrip-20200310');
                ft_defaults
                %         allMicro = '/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/data/microLFP/';
                allMicro = '/castles/nr/projects/h/hanslmas-ieeg-compute\Luca\data\microLFP_dmeanOrth\'; % demeaned and orth. LFP
        
                load('/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/data/allSbj/allSpksHZ.mat', 'allSpks');
end

prevLFPname = [];
for su = 1 : size(allSpks,2)
    disp(su)
    
    allSpks(su).encPhs = [];
    allSpks(su).retPhs = [];
    
    spks = allSpks(su).spks;
    spks = round(spks);
    
    %     if size(spks,1) > 1000
    %         idx = randperm(size(spks, 1))
    %         spks = spks(idx);
    %         spks = spks(1:1000);
    %     end
    
    numSpks    = size(spks,1);
    iu         = allSpks(su).iu;
    bidsID     = allSpks(su).bidsID;
    sesh       = allSpks(su).sesh;
    wirename   = allSpks(su).wirename;
    bundlename = allSpks(su).bundlename;
    
    lfpName = [allMicro, bidsID, '_', sesh, '_', bundlename, '_onlyMicroLFP_dmeanOrth_1000DS_noSPKINT.mat'];
    
    %% ONLY LOAD IN LFP IF IT IS A NEW ONE
    if ~strcmp(lfpName, prevLFPname)
        loaded = 0;
        while loaded == 0
            try
                load(lfpName, 'phsAn');
                loaded = 1;
            catch
                pause(60)
            end
        end
        
        prevLFPname = lfpName;
    end
    
%     %% SELECT BUNDLE LFP
%     selChan      = contains(lfp.label, bundlename);
%     cfg          = [];
%     cfg.channel  = lfp.label(selChan);       % all 8 channels
%     bundLFP      = ft_selectdata(cfg, lfp);  % select
%     
%     numChan      = size(bundLFP.label,1);
    
    %% WHICH TRIALS DOES THIS SU INDEX?
    idxTrl = allSpks(su).idxTrl;
    
    %% IS THIS A (p)IU?
    isIU = allSpks(su).iu>0;
    
    %% ENCODING AND RETRIEVAL TRIGGER
    encTrig  = round(allSpks(su).encTrigger(allSpks(su).hitsIdx,[1 3])*1000);
    retTrig  = round(allSpks(su).retTrigger(allSpks(su).hitsIdx,[1 3])*1000);
    
%     %% REDEFINE TRIALS ACCORDING TO encTrig
%     cfg            = [ ];
%     cfg.trl        = [encTrig(:,1) encTrig(:,2) zeros(size(encTrig,1),1)];
%     bundPHSencTrl  = ft_redefinetrial(cfg, phsAn);
%     
%     %% REDEFINE TRIALS ACCORDING TO retTrig
%     cfg            = [ ];
%     cfg.trl        = [retTrig(:,1) retTrig(:,2) zeros(size(retTrig,1),1)];
%     bundPHSretTrl  = ft_redefinetrial(cfg, phsAn);
    
    %% LOOP OVER TRIALS
    for trl = 1 : size(encTrig,1)        
        % SPIKES THAT OCCUR WITHIN THIS TRIAL DURING ENCODING
        encSpkTrl = spks(spks>=encTrig(trl, 1) & spks <= encTrig(trl, 2)); % ALL SPIKES THAT OCCUR WITHIN THE CURRENT TRIAL
%         if ~isempty(encSpkTrl); encSpkTrl = encSpkTrl(1); end % take the
%         first spike (outdated)
        
        % SPIKES THAT OCCUR WITHIN THIS TRIAL DURING RETRIEVAL
        retSpkTrl = spks(spks>=retTrig(trl, 1) & spks <= retTrig(trl, 2));
%         if ~isempty(retSpkTrl); retSpkTrl = retSpkTrl(1); end % take the
%         first spike (outdated)
        
        allSpks(su).encPhs(trl,1) = {phsAn.phs(:,:,encSpkTrl)}; % in a cell because variable number of spikes per encoding trial
        allSpks(su).retPhs(trl,1) = {phsAn.phs(:,:,retSpkTrl)};
    end
    
    %% WHICH CHANNEL IS THE CHOSEN ONE?    
    chanPhsEnc = cat(3, allSpks(su).encPhs{:});
    chanPhsRet = cat(3, allSpks(su).retPhs{:});

    for chan = 1:8
        for freq = 1:size(chanPhsEnc,2)
            [~, zlvlEnc(chan,freq)] = circ_rtest(squeeze(chanPhsEnc(chan, freq,:))); 
            [~, zlvlRet(chan,freq)] = circ_rtest(squeeze(chanPhsRet(chan, freq,:)));
        end
    end
    
%     for chan = 1:8
%         [plvl(chan,1) zlvl(chan,1)] = circ_rtest(chanPhsEnc);
%         [plvl(chan,2) zlvl(chan,2)] = circ_rtest(chanPhsRet);
%     end
    
    encRetMax = max(zlvlEnc, zlvlRet); % take maximum z-value of encoding vs. retrieval
    chanMax = max(encRetMax, [], 2); % for each channel, take the maximum value over the frequency spectrum
    [~, chosChan] =max(chanMax);
    allSpks(su).LFPchan = chosChan;

end

save('X:\Luca\data\allSbj\allSpksHZ_phs', 'allSpks');