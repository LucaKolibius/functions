%% THIS FUNCTION TAKES THE MICRO-LFP (DS to 1000HZ, BUT NO SPK-INT) AND THEN (1) DEMEANS IT AND THEN (2) ORTHOGONALIZES IT
%  IT ALSO GETS KILLED ON THE VM
function dmeanOrthLFP(server)

switch server
    case 0 %% LAB PC
        
        load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks');
        savefolder = 'X:\Luca\data\microLFP_dmeanOrth\';
        
        
    case 1 %% CASTLES VM
        cd('/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/functions');
        addpath('/castles/nr/projects/h/hanslmas-ieeg-compute/Common/toolboxes/fieldtrip-20200310');
        ft_defaults
        
        load('/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/data/allSbj/allSpksHZ.mat', 'allSpks');
        savefolder = '/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/data/microLFP_dmeanOrth/';
end

for spk = 1 :  size(allSpks,2)
    disp(spk);
    
    %% SKIP REPEATING MICROWIRE BUNDLES
    if spk > 1
        if and(contains(allSpks(spk-1).bidsID, allSpks(spk).bidsID), strcmp(allSpks(spk-1).sesh, allSpks(spk).sesh) ) % same subject + session
            continue
        end
    end
    
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
    curBund = allSpks(spk).bundlename;
    allBund = {allSpks.bundlename};
    
    %% LOAD LFP
    switch server
        case 0 % MY LAPTOP
            load(['X:/Luca/data/microLFP/', bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
            
        case 1 % SERVER
            load(['/castles/nr/projects/h/hanslmas-ieeg-compute/Luca/data/microLFP/', bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    end
    
    %% DEMEAN
    cfg        = [];
    cfg.demean = 'yes';
    data_dm       = ft_preprocessing(cfg, data);
    
    %% ORTHOGONALIZE
    [data_dmOrth, ~, allBund] = orthogonVec(data_dm); % with orthogonalization
    
    %% LOOP OVER BUNDLES
    for bund = 1 : size(allBund,1)
        
        % INDEX LFP OF CURRENT BUNDLE
        curBund = cellfun(@(x) regexp(x,allBund(bund)), {data_dmOrth.label}, 'un', 0);
        curBund = curBund{1};
        curBund = ~cellfun(@isempty, curBund);
        curBund = find(curBund == 1);
        
        % SELECT ALL MW LFP FROM CURRENT BUNDLE
        cfg = [];
        cfg.channel      = data_dmOrth.label(curBund);
        data_dmOrth_bund = ft_selectdata(cfg, data_dmOrth);
        
        %% EXTRACT FREQSPEC
        cfgtf        = [];
        cfgtf.method = 'wavelet';
        cfgtf.width  = 8;
        cfgtf.toi    = 'all';
        cfgtf.foi    = logspace(log10(1),log10(140));
        cfgtf.output = 'fourier';
        fspec        = ft_freqanalysis(cfgtf,data_dmOrth_bund);
        
        phsAn.phs    = squeeze(angle(fspec.fourierspctrm)); % trl(=1) x chan x freq x time
        phsAn.label  = fspec.label;
        phsAn.time   = fspec.time;
        phsAn.freq   = fspec.freq;
        
        clear fspec
        
        %% SAVE LFP
        save([savefolder, bidsID, '_', regexprep(sesh, 'S1b', 'S1'),'_', allBund{bund}, '_onlyMicroLFP_dmeanOrth_1000DS_noSPKINT.mat'], 'data_dmOrth_bund', 'phsAn', '-v7.3');
        
    end % END OF BUNDLE LOOP
end % END OF SPK FOR LOOP
end % END OF FUNCTION
