function [p,num] = findIndexSU_hz(plvl, dynTHmin, encRetTH, afterCue, afterResp)

% encRet min @2 and plvl at 95 seemed to have worked nicely
% maybe dynTH at 1.96 and encRet @1 was sign @6IU
try
    addpath(genpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\functions'))
    load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp.mat', 'allSpks')
catch
    load('/analyse/Project0309/Luca/data/allSbj/allSpksHZ_encTwo_retResp.mat', 'allSpks')
    addpath(genpath('/analyse/Project0309/Luca/functions'))
end

nperm = 10000;
ewpPerm = cell(size(allSpks));

% plvl = 99;
% dynTHmin = 1.96;
% encRetTH = 1;
% % afterCue = 3000; % 3s after cue
% afterCue = 0; % from cue onwards
% afterResp = 0;
allSpks = rmfield(allSpks, 'iu');
allSpks = rmfield(allSpks, 'idxTrl');

for spk = 1 : length(allSpks) 
  
    subj = sub_ID_conversion(allSpks(spk).bidsID, 'yes');
    sesh = allSpks(spk).sesh;
    
    if strcmp(subj, 'P7_ERL')
        subj = 'P07ERL';
    end
    
    try
        cd(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data', filesep, subj, filesep, sesh])
    catch
        cd(['/analyse/Project0309/Luca/data', filesep, subj, filesep, sesh])
    end
    
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
    
    %     %% OLD: WHOLE TRIAL
    %     encTrigger = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1 3])*1000);
    %     retTrigger = round(allSpks(spk).retTrigger(allSpks(spk).hitsIdx,[1 3])*1000);
    
    % ENC: CUE-RESP || RET: CUE-RESP
    encTrigger = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1 3])*1000);     % STIMULUS ONSET UNTIL RESPONSE
    retTrigger = round(allSpks(spk).retTrigger(allSpks(spk).hitsIdx,[1])*1000);       % CUE ONSET UNTIL...
    retTrigger = [ retTrigger round(allSpks(spk).retRT(allSpks(spk).hitsIdx)*1000) ]; % ... RESPONSE
    
    % EXCLUDE CN PERIOD
    encTrigger(:,1) = encTrigger(:,1) + 2000;
%     retTrigger(:,1) = retTrigger(:,1) + 500;
    
    spkTrlEnc = [];
    spkTrlRet = [];
    for trl = 1:length(allSpks(spk).hitsIdx)
        %% ENCODING
        trlLen = (encTrigger(trl,2) - encTrigger(trl,1)) / 1000;
        spkTrl = spkTms(spkTms>=encTrigger(trl,1)+afterCue & spkTms<encTrigger(trl,2)+1000*afterResp) - encTrigger(trl,1) + 1;
        spkTrl = size(spkTrl,1);
        spkTrlEnc(trl) = spkTrl / trlLen;
        
        %% RETRIEVAL
        trlLen = (retTrigger(trl,2) - retTrigger(trl,1))  / 1000;
        spkTrl = spkTms(spkTms>=retTrigger(trl,1)+afterCue & spkTms<retTrigger(trl,2)+1000*afterResp) - retTrigger(trl,1) + 1;
        spkTrl = size(spkTrl,1);
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
        encRetMin  = and(curEncPerm >= encRetTH, curRetPerm >= encRetTH);
        
        % EWP MINIMUM
        ewpMin = curPerm >= dynTH(spk);
        
        % high enough ewq AND high enough firing in encoding and retrieval
                minFiring = and(encRetMin, ewpMin);
%         minFiring = ewpMin;
        
        if any(minFiring)
            
            % TEST FOR INDEXING
            if any(minFiring(ffIdx{spk}) >= 1) & any(minFiring(ppIdx{spk}) >= 1)
                
                permIU(perm) = permIU(perm) + 1;
                permGU(perm) = permGU(perm) + 1;
                
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
% resIU.p = 1 - sum(resIU.num >= permIU) / nperm;
resIU.p = mean(permIU > resIU.num)

resGU.num = sum([allSpks.iu] == 1 | [allSpks.iu] == 2);
resGU.higherIter = sum(resGU.num >= permGU);
% resGU.p = 1 - sum(resGU.num >= permGU) / nperm;
resGU.p = mean(permGU > resGU.num)

inputVar.plvl      = plvl;
inputVar.dynTHmin  = dynTHmin;
inputVar.encRetTH  = encRetTH;
inputVar.afterCue  = afterCue;
inputVar.afterResp = afterResp;

try
    save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp_noCN.mat', 'allSpks', 'resIU', 'resGU', 'inputVar', 'permIU', 'permGU' );
    cd('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj')
catch
    save('/analyse/Project0309/Luca/data/allSbj/allSpksHZ_encTwo_retResp_noCN.mat', 'allSpks', 'resIU', 'resGU', 'inputVar', 'permIU', 'permGU' );
    cd('/analyse/Project0309/Luca/data/allSbj')
end

disp(resIU.p);

p = resIU.p;
num = resIU.num;
end % END OF FUNCTION