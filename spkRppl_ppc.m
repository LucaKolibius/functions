%% SPIKE PPC FOR RIPPLE FREQUENCY BAND (spike ripple ppc)
%  Now adapted to run in a parfor loop
%  This does not take any artefacts into account (they have to be NaN'ed before 

server = 0;

switch server
    % RUN ON MY OWN LAPTOP
    case 0 
        addpath('X:\Common\toolboxes\fieldtrip-20200310');
        ft_defaults
        allMicro = 'X:\Luca\data\microLFP\';
        load('X:\Luca\data\allSbj\allSpks.mat', 'allSpks');
        
    % RUN ON THE SERVER
    case 1
        addpath('/castles/nr/projects/h/hanslmas-ieeg-compute/Common/toolboxes/fieldtrip-20200310');
        ft_defaults
        allMicro = '/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/data/microLFP/';
        load('/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/fileexchange covid/allSpks.mat', 'allSpks');
end

parfor it = 1 : size(allSpks,2)
    tic
    
    spks    = allSpks(it).spks;
    spks = round(spks);
    
    if size(spks,1) > 1000
        idx = randperm(size(spks, 1))
        spks = spks(idx);
        spks = spks(1:1000);
    end
    
    numSpks    = size(spks,1);
    iu         = allSpks(it).iu;
    bidsID     = allSpks(it).bidsID;
    sesh       = allSpks(it).sesh;
    wirename   = allSpks(it).wirename;
    bundlename = allSpks(it).bundlename;
    
    % LOAD LFP
    data = load([allMicro, bidsID, '_', sesh, '_onlyMicroLFP_RAW_1000DS_SPKINT.mat'], 'data');
    data = data.data;

    % SELECT BUNDLE LFP
    selChan      = contains(data.label, bundlename);
    cfg          = [];
    cfg.channel  = data.label(selChan);       % all 8 channels
    data         = ft_selectdata(cfg, data);  % select
    
    numChan      = size(data.label,1);

    % FREQUENCY ANALYSIS TO GET THE PHASE VALUES
    cfgtf        = [];
    cfgtf.method = 'wavelet';
    cfgtf.width  = 8;
    cfgtf.toi    = 1:100; %'all';
    cfgtf.foi    = logspace(log10(1),log10(140))
    cfgtf.output = 'fourier';
    fspec        = ft_freqanalysis(cfgtf,data);
    
    phs          = fspec.fourierspctrm;
    phs          = angle(phs);   % trl=1 / chan / foi / toi
    phs          = squeeze(phs); %       chan / foi / toi
    
    
    % LOOP OVER CHANNELS
    chanFreqDum = [];
    for chan = 1 : numChan
        % ANGLE VALUES FOR CURRENT CHANNEL
        chanPhs = squeeze(phs(chan,:,:)); % foi / toi
        
        freqDum = [];
        
        % LOOP OVER FREQUENCIES
        for curFreq = 1 : size(cfgtf.foi,2)
            
            % ANGLE VALUES FOR CURRENT CHANNEL AND FREQUENCY
            chanFoiPhs = squeeze(chanPhs(curFreq,:)); % toi
            spkDum = [];
            
            % LOOP OVER SPIKES
            for curSpk = 1 : size(spks,1) - 1
                spkDum(curSpk) = nansum( cos( chanFoiPhs ( spks(curSpk)) - chanFoiPhs ( spks(curSpk+1 : end)) ) );
            end
            freqDum(curFreq) = nansum(spkDum) / (numSpks*(numSpks-1)/2);
            
        end
        chanFreqDum(:,chan) = freqDum;
    end
    
    allSpks(it).ppc = chanFreqDum; %% SAVE OUTPUT IN MY STRUCTURE
    
    doneSecs = toc;
    fprintf('Spike number %d. Took %.2f seconds.\n', it, doneSecs)
    fprintf('Spike number %d. Took %.2f seconds.\n', it, doneSecs)
    fprintf('Spike number %d. Took %.2f seconds.\n', it, doneSecs)
    fprintf('Spike number %d. Took %.2f seconds.\n', it, doneSecs)
  
end
save('spkRpplPPC', 'allSpks');