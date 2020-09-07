function [p] = findIndexSU_hz(plvl, encRetTH, afterCue, afterResp)

% encRet min @2 and plvl at 95 seemed to have worked nicely
% maybe dynTH at 1.96 and encRet @1 was sign @6IU

load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks')
nperm = 10000;
ewpPerm = cell(size(allSpks));

% plvl = 99;
% dynTHmin = 1.96;
% encRetTH = 1;
% % afterCue = 3000; % 3s after cue
% afterCue = 0; % from cue onwards
% afterResp = 0;

for spk = 1 : length(allSpks)
    
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
        spkTrl = spkTms(spkTms>=encTrigger(trl,1)+afterCue & spkTms<encTrigger(trl,2)+1000*afterResp) - encTrigger(trl,1) + 1;
        spkTrl = sum(spkTrl);
        spkTrlEnc(trl) = spkTrl / trlLen;
        
        %% RETRIEVAL
        trlLen = (retTrigger(trl,2) - retTrigger(trl,1))  / 1000;
        spkTrl = spkTms(spkTms>=retTrigger(trl,1)+afterCue & spkTms<retTrigger(trl,2)+1000*afterResp) - retTrigger(trl,1) + 1;
        spkTrl = sum(spkTrl);
        spkTrlRet(trl) = spkTrl / trlLen;
    end
    
    %% standardize
    spkTrlEnc = (spkTrlEnc - mean(spkTrlEnc)) / std(spkTrlEnc);
    spkTrlRet = (spkTrlRet - mean(spkTrlRet)) / std(spkTrlRet);
    
    ewp = spkTrlEnc .* spkTrlRet;
    
    %% RESET PREVIOUS ANALYSIS
    allSpks(spk).iu = 0;
    allSpks(spk).idxTrl = false(size(ewp));
    
    %% FIRST LEVEL PERMUTATION
    ewpPerm{spk} = zeros(length(spkTrlEnc), nperm);
    for perm = 1 : nperm
        encPerm{spk}(:,perm) = spkTrlEnc(randperm(size(spkTrlEnc,2)));
        retPerm{spk}(:,perm) = spkTrlRet(randperm(size(spkTrlRet,2)));
        
        ewpPerm{spk}(:,perm) = encPerm{spk}(:,perm) .* retPerm{spk}(:,perm);
    end
    
    dynTH(spk) = prctile(ewpPerm{spk}(:), plvl);
    
    %% fixTH
    if dynTH(spk) < dynTHmin
        dynTH(spk) = dynTHmin;
    end
    
    %% THE TRIAL IS EITHER A FF TRIAL OR A PP TRIAL
    %  ffppTRL = or(ffIdx{spk}, ppIdx{spk});
    
    %% enc+ret minimum
    encRetMin = and(spkTrlEnc >= encRetTH, spkTrlRet >= encRetTH);
    
    %% ewp test
    ewpMin = ewp >= dynTH(spk);
    
    %% high enough firing during encoding, retrieval and ewp
    minFiring = and(encRetMin, ewpMin);
%     minFiring = ewpMin;
    
    if any(minFiring)
        
        % TEST FOR INDEXING
        if any(minFiring(ffIdx{spk}) >= 1) & any(minFiring(ppIdx{spk}) >= 1)
            
            allSpks(spk).iu = 2;
            allSpks(spk).idxTrl = minFiring;
            
        % TEST FOR PUTATIVE INDEXING
        else
            
            allSpks(spk).iu = 1;
            allSpks(spk).idxTrl = minFiring;
            
        end
        
    end
    
end % END OF FIRST LEVEL PERMUTATION

% SECOND LEVEL PERMUTATION
permIU = zeros(1, nperm);
permGU = zeros(1, nperm);
for perm = 1 : nperm
    disp(perm)
    tic
    for spk = 1 : length(allSpks)

        thisPerm = randperm(nperm,1); % chose one out of all nperm (e.g. 10.000) permutations)
        curPerm  = ewpPerm{spk}(:,thisPerm); % this is the current permutation (ewp)
        
        % ENC RET MINIMUM
        curEncPerm = encPerm{spk}(:,thisPerm);
        curRetPerm = retPerm{spk}(:,thisPerm);
        encRetMin = and(curEncPerm >= encRetTH, curRetPerm >= encRetTH);
        
        % EWP MINIMUM
        ewpMin = curPerm >= dynTH(spk);
        
        % high enough ewq AND high enough firing in encoding and retrieval
                minFiring = and(encRetMin, ewpMin);
%         minFiring = ewpMin;
        
        if any(minFiring)
            
            % TEST FOR INDEXING
            if any(minFiring(ffIdx{spk}) >= 1) & any(minFiring(ppIdx{spk}) >= 1)
                
                permIU(perm) = permIU(perm) + 1;
                
            % TEST FOR NUMBER OF GU UNDER THE NULL
            else
                permGU(perm) = permGU(perm) + 1;
            end
            
        end
        
    end % END OF SPIKES
    toc   
end % END OF PERMUTATION

resIU.num = sum([allSpks.iu] == 2);
resIU.higherIter = sum(resIU.num >= permIU);
resIU.p = 1 - sum(resIU.num >= permIU) / nperm;

resGU.num = sum([allSpks.iu] == 1);
resGU.higherIter = sum(resGU.num >= permGU);
resGU.p = 1 - sum(resGU.num >= permGU) / nperm;

save('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks', 'resIU', 'resGU');
cd('X:\Luca\data\allSbj\')
disp(resIU.p);

p = resIU.p;
end % END OF FUNCTION