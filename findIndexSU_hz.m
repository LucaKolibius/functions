function findIndexSU_hz

load('X:\Luca\data\allSbj\allSpks.mat', 'allSpks')
nsecs = 0;
nperm = 10000;
ewpPerm = cell(size(allSpks));
plvl = 99;

for spk = 1 : length(allSpks)
    if isempty(allSpks(spk).spks)
        continue
    end
    
    subj = sub_ID_conversion(allSpks(spk).bidsID, 'yes');
    sesh = allSpks(spk).sesh;
    
    if strcmp(subj, 'P7_ERL')
        subj = 'P07ERL';
    end
    
    cd(['X:/Luca/data', filesep, subj, filesep, sesh])
    abc = dir; cd(abc(3).name);
    p2d = [cd, filesep];
    
    [tableTemplate, ~, ~, ~, ~, ~, ~, ~, ~] = loadLogs(p2d, 1);
    ffIdx{spk}   = ~cellfun(@isempty, strfind(tableTemplate{1,:}, 'ff'));
    ppIdx{spk}   = ~cellfun(@isempty, strfind(tableTemplate{1,:}, 'pp'));
    fpIdx{spk}   = ~cellfun(@isempty, strfind(tableTemplate{1,:}, 'fp'));
    missIdx      =  cellfun(@isempty, strfind(tableTemplate{2,:}, 'hit')); % if hitMiss is 'hit' himiIdx will have a 1 if a trial was a miss
    
    ffIdx{spk}(missIdx) = [];
    ppIdx{spk}(missIdx) = [];
    fpIdx{spk}(missIdx) = [];
    
    %% ENCODING
    spkTms = allSpks(spk).spks;
    
    encTrigger = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1 3])*1000);
    retTrigger = round(allSpks(spk).retTrigger(allSpks(spk).hitsIdx,[1 3])*1000);
    
    spkTrlEnc = [];
    spkTrlRet = [];
    for trl = 1:length(allSpks(spk).hitsIdx)
        %% ENCODING
        trlLen = (encTrigger(trl,2) - encTrigger(trl,1)) / 1000;
        spkTrl = spkTms(spkTms>=encTrigger(trl,1) & spkTms<encTrigger(trl,2)+1000*nsecs) - encTrigger(trl,1) + 1;
        spkTrl = sum(spkTrl);
        spkTrlEnc(trl) = spkTrl / trlLen;
        
        %% RETRIEVAL
        trlLen = (retTrigger(trl,2) - retTrigger(trl,1))  / 1000;
        spkTrl = spkTms(spkTms>=retTrigger(trl,1) & spkTms<retTrigger(trl,2)+1000*nsecs) - retTrigger(trl,1) + 1;
        spkTrl = sum(spkTrl);
        spkTrlRet(trl) = spkTrl / trlLen;
    end
    
    ewp = spkTrlEnc .* spkTrlRet;
    
    %% RESET PREVIOUS ANALYSIS
    allSpks(spk).iu = 0;
    allSpks(spk).idxTrl = false(size(ewp));
    
    %% FIRST LEVEL PERMUTATION
    ewpPerm{spk} = zeros(length(spkTrlEnc), nperm);
    for perm = 1 : nperm
        encPerm = spkTrlEnc(randperm(size(spkTrlEnc,2),1));
        retPerm = spkTrlRet(randperm(size(spkTrlRet,2),1));
        
        ewpPerm{spk}(:,perm) = encPerm .* retPerm;
    end
    
    dynTH(spk) = prctile(ewpPerm{spk}(:), plvl);
    
    % THE TRIAL IS EITHER A FF TRIAL OR A PP TRIAL
    ffppTRL = or(ffIdx{spk}, ppIdx{spk});
    
    % TEST FOR INDEXING
    if any(ewp(ffIdx{spk}) >= dynTH(spk)) & any(ewp(ppIdx{spk}) >= dynTH(spk))
        
        allSpks(spk).iu = 2;
        allSpks(spk).idxTrl = and( ewp >= dynTH(spk), ffppTRL); % FFPP TRIALS THAT ALSO EXCEED THE FIRING THRESHOLD
        
        % TEST FOR PUTATIVE INDEXING
    elseif any(ewp >= dynTH(spk))
        allSpks(spk).iu = 1;
        allSpks(spk).idxTrl = ewp >= dynTH(spk);
    end
    
end

% SECOND LEVEL PERMUTATION
permIU = zeros(1, nperm);
permGU = zeros(1, nperm);
for perm = 1 : nperm
    for spk = 1 : length(allSpks)
        
        if isempty(allSpks(spk).spks)
            continue
        end
        
        thisPerm = randperm(nperm,1);
        curPerm  = ewpPerm{spk}(:,thisPerm);
        
        % TEST FOR NUMBER OF IU UNDER THE NULL
        if any(curPerm(ffIdx{spk}) >= dynTH(spk) & any(curPerm(ppIdx{spk}) >= dynTH(spk)))
            permIU(perm) = permIU(perm) + 1;
        end
        
        % TEST FOR NUMBER OF GU UNDER THE NULL
        if any(curPerm >= dynTH(spk))
            permGU(perm) = permGU(perm) + 1;
        end
        
    end
end
resIU.num = sum([allSpks.iu] == 2);
resIU.higherIter = sum(resIU.num > permIU);
resIU.p = 1 - sum(resIU.num > permIU) / nperm;

resGU.num = sum([allSpks.iu] == 1);
resGU.higherIter = sum(resGU.num > permGU);
resGU.p = 1 - sum(resGU.num > permGU) / nperm;

save('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks', 'resIU', 'resGU');
end % END OF FUNCTION