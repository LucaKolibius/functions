load('//analyse4.psy.gla.ac.uk/project0309/Luca/data/allSbj/allSpksHZ_encTwo_retResp_noCN.mat', 'allSpks');
for su = 1:length(allSpks)
    if length(allSpks(su).spks) < 1000
        continue
    end
    
    bidsID = allSpks(su).bidsID;
    sesh   = allSpks(su).sesh;
    chan   = allSpks(su).wirename;
    
    abc = dir(['\\analyse4.psy.gla.ac.uk\project0309\Data Continuous\', bidsID, '_', sesh, '_micro', '*']);
    load([abc.folder, filesep, abc.name], 'data_micro', 'spike_waveform')
    
    cfg = [];
    cfg.resamplefs = 1000;
    cfg.demean     = 'yes';
    cfg.detrend    = 'yes';
    data = ft_resampledata(cfg, data_micro);
    
    
    ypar(1) = min(data.trial{1}(:,1:10000),[], 'all');
    ypar(2) = max(data.trial{1}(:,1:10000),[], 'all');
    
    cfg = [];
    cfg.continuous = 'yes';
    cfg.viewmode = 'vertical';
    cfg.channel = chan;
    % cfg.ylim = [ypar(1) ypar(2)]
    cfg.ylim = [-50 50];
    tt  = ft_databrowser(cfg, data);
    
end

