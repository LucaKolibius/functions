%% THIS FUNCTION TAKES THE MICRO-LFP (DS to 1000HZ, BUT NO SPK-INT) AND THEN (1) DEMEANS IT AND THEN (2) ORTHOGONALIZES IT

function dmeanOrthLFP

load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks');
lfpDir     = dir('X:\Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_noSPKINT.mat');  % NO SPK-INT
savefolder = 'X:\Luca\data\microLFP_dmeanOrth\';

for spk = 1 :  size(allSpks,2)
    
    
    %% SKIP REPEATING MICROWIRE BUNDLES
    if spk > 1
        if and(and(contains(allSpks(spk-1).bidsID, allSpks(spk).bidsID), strcmp(allSpks(spk-1).sesh, allSpks(spk).sesh) ), strcmp(allSpks(spk).bundlename, allSpks(spk-1).bundlename)); % same subject + session
            continue
        end
    end
    
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
    curBund = allSpks(spk).bundlename;
    allBund = {allSpks.bundlename};
    
    %% LOAD LFP
    load([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    
    %% DEMEAN
    cfg        = [];
    cfg.demean = 'yes';
    data       = ft_preprocessing(cfg, data);
    
    %% ORTHOGONALIZE
    data = orthogonVec(data); % with orthogonalization
    
    %% EXTRACT FREQSPEC
    cfgtf        = [];
    cfgtf.method = 'wavelet';
    cfgtf.width  = 8;
    cfgtf.toi    = 'all';
    cfgtf.foi    = logspace(log10(1),log10(140));
    cfgtf.output = 'fourier';
    fspec        = ft_freqanalysis(cfgtf,data);   
    
    %% SAVE LFP
    save([savefolder, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_dmeanOrth_1000DS_noSPKINT.mat'], 'data', 'fspec');
    
end % END OF SPK FOR LOOP
end % END OF FUNCTION
