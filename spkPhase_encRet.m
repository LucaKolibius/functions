%% SPIKE PPC 
%  This does not take any artefacts into account (they have to be NaN'ed before
clear
server = 1;

switch server
    % RUN ON MY OWN LAPTOP
    case 0
        addpath('X:\Common\toolboxes\fieldtrip-20200310');
        addpath('X:\Luca\toolboxes\CircStat2012a');
        ft_defaults
        %         allMicro = 'X:\Luca\data\microLFP\'; % non interpolated LFP
        allMicro = 'X:\Luca\data\microLFP_dmeanOrth\'; % demeaned and orth. LFP
        load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks');
        
        % RUN ON THE SERVER
        %     case 1
        %         addpath('/castles/nr/projects/h/hanslmas-ieeg-compute/Common/toolboxes/fieldtrip-20200310');
        %         ft_defaults
        %         %         allMicro = '/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/data/microLFP/';
        %         allMicro = '/castles/nr/projects/h/hanslmas-ieeg-compute\Luca\data\microLFP_dmeanOrth\'; % demeaned and orth. LFP
        
        %         load('/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/data/allSbj/allSpksHZ.mat', 'allSpks');
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
    
    lfpName = [allMicro, bidsID, '_', sesh, '_onlyMicroLFP_dmeanOrth_1000DS_noSPKINT.mat'];
    
    %% ONLY LOAD IN LFP IF IT IS A NEW ONE
    if ~strcmp(lfpName, prevLFPname)
        lfp = load(lfpName, 'phsAn');
        lfp = lfp.data;
        prevLFPname = lfpName;
    end
    
    %%% SELECT BUNDLE LFP
    selChan      = contains(lfp.label, bundlename);
    cfg          = [];
    cfg.channel  = lfp.label(selChan);       % all 8 channels
    bundLFP      = ft_selectdata(cfg, lfp);  % select
    
    numChan      = size(bundLFP.label,1);
    
    %% WHICH TRIALS DOES THIS SU INDEX?
    idxTrl = allSpks(su).idxTrl;
    
    %% IS THIS A (p)IU?
    isIU = allSpks(su).iu>0;
    
    %% ENCODING AND RETRIEVAL TRIGGER
    encTrig  = round(allSpks(su).encTrigger(allSpks(su).hitsIdx,[1 3])*1000);
    retTrig  = round(allSpks(su).retTrigger(allSpks(su).hitsIdx,[1 3])*1000);
    
    %% REDEFINE TRIALS ACCORDING TO encTrig
    cfg            = [ ];
    cfg.trl        = [encTrig(:,1) encTrig(:,2) zeros(size(encTrig,1),1)];
    bundLFPencTrl  = ft_redefinetrial(cfg, bundLFP);
    
    %% REDEFINE TRIALS ACCORDING TO retTrig
    cfg            = [ ];
    cfg.trl        = [retTrig(:,1) retTrig(:,2) zeros(size(retTrig,1),1)];
    bundLFPretTrl  = ft_redefinetrial(cfg, bundLFP);
    
    %% LOOP OVER TRIALS
    for trl = 1 : size(encTrig,1)        
        % SPIKES THAT OCCUR WITHIN THIS TRIAL DURING ENCODING
        encSpkTrl = spks>=encTrig(1,trl) && spks <= encTrig(2,trl);
        if ~isempty(encSpkTrl); encSpkTrl = encSpkTrl(1); end % take the first spike
        
        % SPIKES THAT OCCUR WITHIN THIS TRIAL DURING RETRIEVAL
        retSpkTrl = spks>=retTrig(1,trl) && spks <= retTrig(2,trl);
        if ~isempty(retSpkTrl); retSpkTrl = retSpkTrl(1); end % take the first spike
        
        allSpks(su).encPhs(:,:,trl) = phsAn.phs(encSpkTrl:::);
        allSpks(su).retPhs(:,:,trl) = phsAn.phs(retSpkTrl:::);
    end
    
    %% WHICH CHANNEL IS THE CHOSEN ONE?    
    chanPhsEnc = resize(allSpks(su.encPhs), 8, [], 1); % CHECK IF CORRECT
    chanPhsRet = resize(allSpks(su.retPhs), 8, [], 1); % CHECK IF CORRECT

    for chan = 1:8
        [plvl(chan,1) zlvl(chan,1)] = circ_rtest(chanPhsEnc);
        [plvl(chan,2) zlvl(chan,2)] = circ_rtest(chanPhsRet);
    end
    
    [encRetVal,~] = min(plvl,[],2); % take either minimum of encoding or retrieval
    [~, chanChos] = min(encRetVal); % take lowest channel
    
    allSpks(su).LFPchan = chanChos;

end

save('X:\Luca\data\allSbj\allSpksHZ_phs', 'allSpks');