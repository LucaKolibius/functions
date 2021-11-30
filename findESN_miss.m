% Main function to detect miss and half miss ESNs. Excludes putative concept neurons.
function findESN_miss

% poolobj  = gcp('nocreate');
% if isempty(poolobj) % If there is no parallel pool, create one
%     parpool('local', 20)
% end

plvl      = 99;
encRetTH  = 1.645;
nperm     = 10000;

folderpath = '\\analyse4.psy.gla.ac.uk\project0309\Luca\ESN_code_data\code_ESN\data'; % wherever you have downloaded this folder
addpath(genpath(folderpath))
load([folderpath, '\allSpks.mat']);

tic
for himi = 2:3
    
    ewpPerm    = cell(size(allSpks));
    isCN       = []; dynTH = []; thCN = [];
    allEncShuf = {}; encPerm = {}; retPerm = {}; spksTr = {};
    
    for spk = 1 : size(allSpks,2)
        timePassed = toc;
        lineLength = fprintf('Analysing SU#%d of %d | Time passed: %.1f min \n', spk, size(allSpks,2), (timePassed/60));
        
        % RESET PREVIOUS ANALYSIS
        switch himi
            
            case 2 % half hits
                
                allSpks(spk).hmESN       = 0;
                allSpks(spk).hmReinstTrl = false(1, sum(allSpks(spk).himiDx == 2));
                
            case 3 % miss
                
                allSpks(spk).mESN       = 0;
                allSpks(spk).mReinstTrl = false(1, sum(allSpks(spk).himiDx == 3));
                
        end
        
        % WHICH TRIALS CONTAINED TWO FACES (ff), TWO PLACES (pp), OR A PLACE AND A FACE AS ASSOCIATE IMAGES (fp)
        %% THIS IS ONLY FOR HITS!! I DONT HAVE IT FOR MISSES
        %     ffIdx{spk} = allSpks(spk).ffIdx;
        %     ppIdx{spk} = allSpks(spk).ppIdx;
        %     fpIdx{spk} = allSpks(spk).fpIdx;
        
        % GET THE SPIKES
        spkTms = allSpks(spk).spks;
        
        % ENC: CUE-RESP || RET: CUE-RESP
        missIdx = allSpks(spk).himiDx == himi;
        encTrigger = round(allSpks(spk).encTrigger(missIdx,[1 3])*1000);     % CUE ONSET UNTIL RESPONSE
        retTrigger = round(allSpks(spk).retTrigger(missIdx,[1])*1000);       % CUE ONSET UNTIL...
        retTrigger = [ retTrigger round(allSpks(spk).retRT(missIdx)*1000) ]; % ... RESPONSE
        
        % EXCLUDE CN PERIOD
        encTrigger(:,1) = encTrigger(:,1) + 2000; % NOW STARTING AT RESPONSE ONSET
        
        % SEGMENT THE SPIKETIMES INTO TRIALS
        spkTrlEnc    = [];
        spkTrlRet    = [];
        
        for trl = 1:sum(missIdx)
            % ENCODING
            trlLen         = (encTrigger(trl,2) - encTrigger(trl,1)) / 1000;
            spkTrl         = spkTms(spkTms>=encTrigger(trl,1) & spkTms<encTrigger(trl,2)) - encTrigger(trl,1) + 1;
            spkTrl         = size(spkTrl,1);
            spkTrlEnc(trl) = spkTrl / trlLen; % TRANSFORM IN Hz
            
            % RETRIEVAL
            trlLen         = (retTrigger(trl,2) - retTrigger(trl,1))  / 1000;
            spkTrl         = spkTms(spkTms>=retTrigger(trl,1) & spkTms<retTrigger(trl,2)) - retTrigger(trl,1) + 1;
            spkTrl         = size(spkTrl,1);
            spkTrlRet(trl) = spkTrl / trlLen;
        end
        
        % Z-SCORE SPIKE NUMBER
        spkTrlEnc = (spkTrlEnc - mean(spkTrlEnc)) / std(spkTrlEnc);
        spkTrlRet = (spkTrlRet - mean(spkTrlRet)) / std(spkTrlRet);
        
        ewp = spkTrlEnc .* spkTrlRet; % ELEMENT WISE PRODUCT
        
        % FIRST LEVEL PERMUTATION
        ewpPerm{spk} = zeros(length(spkTrlEnc), nperm);
        for perm = 1 : nperm
            
            % SET SHUFFLE ORDER FOR ENCODING AND RETRIEVAL
            encShuf                   = randperm(size(spkTrlEnc,2));
            retShuf                   = randperm(size(spkTrlRet,2));
            
            allEncShuf{spk}(perm,:)   = encShuf; % save for second level cueCN permutation test
            
            % SHUFFLE SPIKE NUMBER FOR 1ST AND 2ND LEVEL PERMUTATION FOR ESN IDENTIFICATION
            encPerm{spk}(:,perm) = spkTrlEnc(encShuf);
            retPerm{spk}(:,perm) = spkTrlRet(retShuf);
            
            ewpPerm{spk}(:,perm) = encPerm{spk}(:,perm) .* retPerm{spk}(:,perm);
        end
        
        dynTH(spk) = prctile(ewpPerm{spk}(:), plvl);
        
        % DOES THE FIRING EXCEED THE THRESHOLD (1.645) DURING ENCODING AND RETRIEVAL?
        encRetMin = and(spkTrlEnc >= encRetTH, spkTrlRet >= encRetTH);
        
        % DOES THE REINSTATEMENT VALUE (EWP) EXCEED THE COMPUTED THRESHOLD
        ewpMin = ewp >= dynTH(spk);
        
        % DURING WHICH TRIALS DO WE HAVE ABOVE THRESHOLD REINSTATEMENT AND ABOVE
        % THRESHOLD FIRING IN ENCODING + RETRIEVAL
        minFiring = and(encRetMin, ewpMin);
        
        %% CHECK IF SU IS cueCN
        % RETRIEVAL TRIGGER (CUE)
        cueRet = allSpks(spk).retTrigger(missIdx, [1])*1000;
        % ENCODING TRIGGER (CUE)
        cueEnc = allSpks(spk).encTrigger(missIdx, [1])*1000;
        
        % COMPUTE BASELINES
        spksEncBL = insertSpikeNumber(cueEnc, spkTms, [-1000 -300]);
        spksRetBL = insertSpikeNumber(cueRet, spkTms, [-1000 -300]);
        spksBL    = vertcat(spksEncBL, spksRetBL);
        
        % COMPUTE ACTIVITY DURING CN PERIOD
        spksTr{spk}  = insertSpikeNumber(cueEnc, spkTms, [ 300   1000]);
        spksTrEmp    = spksTr{spk}(minFiring,:);
        
        medFir    = median(spksTrEmp,'all'); % median firing rate during reinstated episodes
        
        blMean    = mean(spksBL);
        blStd     = std(spksBL);
        thCN(spk) = blMean + blStd * 5; % threshold firing is at mean(baseline) + 5*STD(Baseline)
        
        if medFir >= 2 && medFir > thCN(spk) % median firing rate has to be at least two
            isCN(spk) = true;
            fprintf(repmat('\b', 1, lineLength));
            continue
        end
        
        %% CODE RESULTS INTO allSpks
        if any(minFiring)
            
            %         % TEST FOR CONSERVATIVE ESN (WITHOUTH CN!)
            %         if any(minFiring(ffIdx{spk}) >= 1) & any(minFiring(ppIdx{spk}) >= 1)
            %
            %             allSpks(spk).ESN       = 2;
            %             allSpks(spk).reinstTrl = minFiring; % index which episodes are reinstated
            %
            %         % TEST FOR CONSERVATIVE ESN
            %         else
            switch himi
                
                case 2 % himi
                    
                    allSpks(spk).hmESN = 1;
                    allSpks(spk).hmReinstTrl = minFiring;
                    
                case 3 % miss
                    
                    allSpks(spk).mESN       = 1;
                    allSpks(spk).mReinstTrl = minFiring;
                    
            end
            %         end
            
        end
        
        fprintf(repmat('\b', 1, lineLength)); % this is nice because it deletes the old progress update
        
    end % END OF FIRST LEVEL PERMUTATION
    
    
    %% SECOND LEVEL PERMUTATION
    % permESNc = zeros(1, nperm);
    permESN = zeros(1, nperm);
    parfor perm = 1 : nperm
        disp(perm)
        for spk = 1 : length(allSpks)
            
            thisPerm = randperm(nperm,1);        % chose one out of all nperm (e.g. 10.000) permutations)
            curPerm  = ewpPerm{spk}(:,thisPerm); % this is the current permutation (ewp)
            
            % ENC RET MINIMUM
            curEncPerm = encPerm{spk}(:,thisPerm);
            curRetPerm = retPerm{spk}(:,thisPerm);
            encRetMin  = and(curEncPerm >= encRetTH, curRetPerm >= encRetTH);
            
            % EWP MINIMUM
            ewpMin = curPerm >= dynTH(spk);
            
            % high enough ewq AND high enough firing in encoding and retrieval
            minFiring = and(encRetMin, ewpMin);
            
            %% IS THIS PERMUTATION A cueCN?
            thisShuf      = allEncShuf{spk}(thisPerm,:);          % apply the same shuffling to the CN trial values as in this permutation
            spksTrPerm    = spksTr{spk}(thisShuf(minFiring),:); % which trials are reinstated in this permutation?
            
            medFirPerm    = median(spksTrPerm,'all'); % what is the median firing rate of these permuted reinstated trials? would we consider them cueCN?
            
            % CONTINUE WITH NEXT PERMUTATION IF THIS IS ONE IS A cueCN
            if medFirPerm >= 2 && medFirPerm > thCN(spk)
                continue
            end
            
            if any(minFiring)
                
                %             % TEST FOR CONSERVATIVE ESN (REINSTATE ONE FF AND ONE PP TRIAL)
                %             if any(minFiring(ffIdx{spk}) >= 1) & any(minFiring(ppIdx{spk}) >= 1)
                %
                %                 permESNc(perm) = permESNc(perm) + 1;
                %                 permESN (perm) = permESN(perm) + 1;
                %
                %             % TEST FOR NUMBER OF ESN UNDER THE NULL
                %             else
                permESN(perm) = permESN(perm) + 1;
                %             end
                
            end
            
        end % END OF SPIKES
    end % END OF PERMUTATION
    
    
    %% RESULTS OF 2ND ORDER PERMUTATION TEST
    % resESNc.num = sum([allSpks.ESN] == 2);
    % resESNc.higherIter = sum(resESNc.num >= permESNc);
    % resESNc.p = mean(permESNc >= resESNc.num)
    
    switch himi
        case 2 % half hits
            
            HMresESN.num        = sum([allSpks.hmESN] == 1 | [allSpks.hmESN] == 2);
            HMresESN.higherIter = sum(HMresESN.num >= permESN);
            HMresESN.p          = mean(permESN >= HMresESN.num)
            HMresESN.perm       = permESN;
            
        case 3 % miss
            
            MresESN.num        = sum([allSpks.mESN] == 1 | [allSpks.mESN] == 2);
            MresESN.higherIter = sum(MresESN.num >= permESN);
            MresESN.p          = mean(permESN >= MresESN.num)
            MresESN.perm       = permESN;
            
    end
    
end % END OF HIMI

cd(folderpath)
save([folderpath, '\allSpks.mat'], 'allSpks', 'resESNc', 'resESN', 'MresESN', 'HMresESN', 'permESNc', 'permESN');

end % END OF FUNCTION