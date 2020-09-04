clear
lfpDir = dir('X:\Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_noSPKINT.mat.mat');  % CHANGED TO NO SPK INT
artFold = 'Z:\hanslmas-ieeg-compute\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS
load('X:\Luca\data\allSbj\allTrig.mat', 'trigALL')
counter = 1;
for ses = 1 : size(lfpDir,1)
    disp(ses)
    
    %% GET: bidsID + sesh
    bidsID = lfpDir(ses).name; sesh = bidsID;
    bidsID = bidsID(1:8);
    sesh   = sesh(10:11);
    subjID = sub_ID_conversion(bidsID, 1);
    
    %% LOAD IN THE DATA
    load([lfpDir(ses).folder, filesep, lfpDir(ses).name], 'data');
    
    %% DETECT RIPPLES
    [~, rpplsWire, bndLab] = calcRppl (data, bidsID, sesh, artFold);
    
    %% MAKE INTO A NEAT STRUCTURE
    for bund = 1 : size(rpplsWire,1)
        clear elecLoc
        cd('Z:\hanslmas-ieeg-compute\Luca\data\')
        cd(subjID); cd(sesh);  %% 1007 sesh 1b??
        abc = dir('2*'); cd(abc.name);
        cd('advancedAnalysis\elecLoc');
        abc = dir('elecLoc*.mat'); load(abc.name);
        bundname = bndLab(bund);
        idx = strcmp(elecLoc(:,1), bundname);
        
        % WHICH TRIALS IN THAT BUNDLE ARE INDEXED (trlIdxBund)
        allBund  = {trigALL.wirename};
        allBund  = cellfun(@(x) x(1:end-1), allBund, 'un', 0);
        sameBund = and(and(contains({trigALL.bidsID}, bidsID), strcmp({trigALL.sesh}, sesh) ), contains(allBund, bndLab(bund))); % same subject + session + bundle
        idxTrl   = any(vertcat(trigALL(sameBund).idxTrl),1); % these are the trials that are indexed in that wire
        
        % GET ENCODING TRIGGER FROM PREVIOUS VARIABLE
        % incomplete for bundles w/o single units
        idxBids = cell2mat(cellfun(@(x) strcmp(x, bidsID), {trigALL.bidsID}, 'un', 0));
        idxSesh = cell2mat(cellfun(@(x) strcmp(x, sesh), {trigALL.sesh}, 'un', 0));
        idxBoth = and(idxBids, idxSesh);
        idxBoth = find(idxBoth == 1);
        
        rpplBund(counter).bidsID     = bidsID;
        rpplBund(counter).sesh       = sesh;
        rpplBund(counter).bundname   = bndLab{bund};
        rpplBund(counter).rppls      = rpplsWire{bund};
        rpplBund(counter).idxTrlBund = idxTrl;
        rpplBund(counter).location   = elecLoc(idx,2);
        
        if isempty(idxBoth)
            counter  = counter + 1;
            continue
        else
            idxBoth = idxBoth(1);
        end
        
        % if the variable does not have a idxTrlBund it means there were no SU on the bundle. if there is no encTrigger it means there were no SU in the session (IIRC)
        rpplBund(counter).encTrigger = trigALL(idxBoth).encTrigger(trigALL(idxBoth).hitsIdx,:);
        counter  = counter + 1;
        
    end
    
end


%% FOR MACRO
clear
addpath(genpath('X:\Luca\functions\'));
addpath('X:\Common\toolboxes\fieldtrip-20200310'); ft_defaults;
lfpDir = dir('X:\Luca\data\macroLFP\sub-*_onlyMacroLFP_RAW_1000DS.mat');  % CHANGED TO NO SPK INT
% artFold = 'X:\George\Analysis\Artefact Rejection\Data Continuous 1000Hz\';    % ARTEFACTS
load('X:\Luca\data\allSbj\allTrig.mat', 'trigALL')
counter = 1;
rpplBund = [];
macro = 1;
ortho = 0;
for ses = 1 : size(lfpDir,1)
    disp(ses)
    
    %% GET: bidsID + sesh
    bidsID = lfpDir(ses).name; sesh = bidsID;
    bidsID = bidsID(1:8);
    sesh   = sesh(10:11);
    subjID = sub_ID_conversion(bidsID, 1);
    
    %% LOAD IN THE DATA
    load([lfpDir(ses).folder, filesep, lfpDir(ses).name], 'data');
    
    % I had to add this here to transform the macro label convention to the
    % micro label convention. because why not have completely different
    % naming conventions for the same bundle in macro and micro
    if macro == 1
        
        if strcmp(bidsID(5), '0') % bham
            data.label = cellfun(@(x) x(end-6:end-1), data.label, 'un',0);
            for lab = 1 : size(data.label,1)
                curLab = data.label(lab);
                newLab = [];
                
                
                if contains(curLab, 'A')
                    newLab = [newLab, 'ant'];
                elseif contains(curLab, 'PHG')
                    newLab = [newLab, 'para'];
                elseif contains(curLab, 'P')
                    newLab = [newLab, 'post'];
                elseif contains(curLab, 'M')
                    newLab = [newLab, 'mid'];
                end
                
                if contains(curLab, 'H')
                    newLab = [newLab, 'Hipp'];
                end
                
                if contains(curLab, 'L')
                    newLab = [newLab, 'L'];
                    
                elseif contains(curLab, 'R')
                    newLab = [newLab, 'R'];
                end
                
                newLab = [newLab, '1']; % calcRppl expects a wire number
                data.label(lab) = {newLab};
            end % end of case Bham Loop over Bundles
            
        elseif strcmp(bidsID(5), '1') % ERLANGEN
            data.label = cellfun(@(x) x(end-1), data.label, 'un',0);
            for lab = 1 : size(data.label,1)
                curLab = data.label(lab);
                newLab = ['M'];
                
                if strcmp(bidsID, 'sub-1003') % MW in sub-1003 are called e.g. "MiA" instead of "MA"
                    newLab = [newLab, 'i'];
                end
                
                newLab = [newLab, curLab{1}(end)];
                newLab = [newLab, '1']; % calcRppl expects a wire number
                data.label(lab) = {newLab};
            end
            
        end % BHAM OR ERLANGEN
    end % CHECK IF IT IS MACRO RECORDING
    
    %% DETECT RIPPLES
    [~, rpplsWire, bndLab] = calcRppl (data, bidsID, sesh, 1, ortho);
    
    %% MAKE INTO A NEAT STRUCTURE
    for bund = 1 : size(rpplsWire,1)
        clear elecLoc
        cd('X:\Luca\data\')
        cd(subjID); cd(sesh);  %% 1007 sesh 1b??
        abc = dir('2*'); cd(abc.name);
        cd('advancedAnalysis\elecLoc');
        abc = dir('elecLoc*.mat'); load(abc.name);
        bundname = bndLab(bund);
        idx = strcmp(elecLoc(:,1), bundname);
        
        % WHICH TRIALS IN THAT BUNDLE ARE INDEXED (trlIdxBund)
        allBund  = {trigALL.wirename};
        allBund  = cellfun(@(x) x(1:end-1), allBund, 'un', 0);
        sameBund = and(and(contains({trigALL.bidsID}, bidsID), strcmp({trigALL.sesh}, sesh) ), contains(allBund, bndLab(bund))); % same subject + session + bundle
        idxTrl   = any(vertcat(trigALL(sameBund).idxTrl),1); % these are the trials that are indexed in that wire
        
        % GET ENCODING TRIGGER FROM PREVIOUS VARIABLE
        % incomplete for bundles w/o single units
        idxBids = cell2mat(cellfun(@(x) strcmp(x, bidsID), {trigALL.bidsID}, 'un', 0));
        idxSesh = cell2mat(cellfun(@(x) strcmp(x, sesh), {trigALL.sesh}, 'un', 0));
        idxBoth = and(idxBids, idxSesh);
        idxBoth = find(idxBoth == 1);
        
        rpplBund(counter).bidsID     = bidsID;
        rpplBund(counter).sesh       = sesh;
        rpplBund(counter).bundname   = bndLab{bund};
        rpplBund(counter).rppls      = rpplsWire{bund};
        rpplBund(counter).idxTrlBund = idxTrl;
        rpplBund(counter).location   = elecLoc(idx,2);
        
        if isempty(idxBoth)
            counter  = counter + 1;
            continue
        else
            idxBoth = idxBoth(1);
        end
        
        % if the variable does not have a idxTrlBund it means there were no SU on the bundle. if there is no encTrigger it means there were no SU in the session (IIRC)
        rpplBund(counter).encTrigger = trigALL(idxBoth).encTrigger(trigALL(idxBoth).hitsIdx,:);
        counter  = counter + 1;
        
    end
    
end
save('X:\Luca\data\allSbj\rpplBund_macro.mat', 'rpplBund');
    

%% ANALYSIS PART
clear
load('X:\Luca\data\allSbj\allSpks.mat', 'allSpks');
load('X:\Luca\data\allSbj\rpplBund_macro.mat', 'rpplBund');
numMW = 1; % for macro recording
numMW = 8; % for micro recording
for tw = 0:1:3
    for maxMW = 1 % 0:1:1
        
        numRpplIdx = [];
        numRpplNdx = [];
        lenRpplNdx = [];
        lenRpplIdx = [];
        rpplDenIdx = [];
        rpplDenNdx = [];
        coincIdx   = [];
        coincNdx   = [];
        
        for bund = 1 : length(rpplBund)
            
            bidsID   = rpplBund(bund).bidsID;
            sesh     = rpplBund(bund).sesh;
            bundname = rpplBund(bund).bundname;
            idxTrl   = rpplBund(bund).idxTrlBund; % INDEXED TRIALS

            % IF WE HAVE NO LOCATION INFORMATION, CONTINUE
            if isempty(rpplBund(bund).location)
                continue
            end
            
            % IF BUNDLE HAS NO SPIKES OR IS NOT IN THE HIPPOCAMPUS, CONTINUE WITH
            % THE NEXT ONE
            if isempty(rpplBund(bund).idxTrlBund) || ~strcmp(rpplBund(bund).location, 'hipp')
                continue
            end
            
            idxSpk = and(and(strcmp({allSpks.bidsID}, bidsID), strcmp({allSpks.sesh}, sesh)), strcmp({allSpks.bundlename}, bundname));
            idxSpk = find(idxSpk == 1);
            
            % DEFINE TIME WINDOW
            encTrig = rpplBund(bund).encTrigger(:, [1, 3]);
            encTrig = round(encTrig * 1000);
            
            % DIFFERENT TIME WINDOWS OF INTEREST
            switch  tw
                case 0                                                           % WHOLE TRIAL
                    encTrig = [encTrig(:,1) encTrig(:,2)];
                case 1                                                           % PRE-CUE
                    encTrig = [encTrig(:,1)-1000 encTrig(:,1)];
                case 2                                                           % PERI-CUE
                    encTrig = [encTrig(:,1)-1000 encTrig(:,1)+1000];
                case 3                                                           % PERI-RESP
                    encTrig = [encTrig(:,2)-1000 encTrig(:,2)+1000];
            end
            
            trlNum       = length(rpplBund(bund).idxTrlBund);
            
            numRppl      = zeros(trlNum,8);
            rpplDen      = zeros(trlNum,8);
            lenRppl      = zeros(trlNum,8);
            rpplSpkCoinc = zeros(trlNum,8);
            
            for mw = 1:numMW
                rppls = rpplBund(bund).rppls{mw}; % START + END OF THE RIPPLES ON THE CURRENT MICROWIRE
                    
                if isempty(rppls) % NOT SURE IF THIS IS OK
                    continue
                end
                
                % LOOP OVER TRIALS
                for trl = 1 : size(encTrig,1)
                    lenRpplAll = [];
                    trlRpplDen = zeros(1, [encTrig(trl,2)-encTrig(trl,1)+1]); % ripple density in that trial
                    

                    % LOOP OVER RIPPLES
                    for ripNum = 1 : size(rppls,1)
                        
                        % IS THERE AN OVERLAP OF THE CURRENT RIPPLE WITH THE
                        % CURRENT TIME WINDOW?
                        if any(intersect([rppls(ripNum,1):rppls(ripNum,2)], [encTrig(trl,1) : encTrig(trl,2)]))                                     % is there any overlap with the current ripple and the current trial?
                          
                            numRppl(trl,mw) = numRppl(trl,mw) +1;                                                                                   % one more ripple occurs within the current trial
                            lenRpplAll = [ lenRpplAll, length(rppls(ripNum,1): rppls(ripNum,2)) ];                                                  % length of the current ripple. is later averaged over all ripples in that trial
                            trlRpplDen(intersect([rppls(ripNum,1):rppls(ripNum,2)], [encTrig(trl,1) : encTrig(trl,2)]) -encTrig(trl,1) +1) = 1;     % ripple density. is computed once all ripples have been considered for that trial.
                        end
                        
                    end % END OF RIPPLE LOOP
                    
                    for spk = 1 : size(idxSpk,2)
                        
                        spkTms    = allSpks(idxSpk(spk)).spks; % all spikes on the current bundle
                        spkTrls   = spkTms(spkTms >= encTrig(trl,1) & spkTms <= encTrig(trl,2)) - encTrig(trl,1)+1; % the spikes that occur in this trial
                        temp(spk) = mean(trlRpplDen(spkTrls)) / mean(trlRpplDen);
                   
                    end % END OF SPIKE LOOP 
                        % outside of ripple loop!
                    
                    lenRppl(trl,mw)      = mean(lenRpplAll);
                    rpplDen(trl,mw)      = mean(trlRpplDen);
                    rpplSpkCoinc(trl,mw) = mean(temp);
                    
                end % END OF TRIAL LOOP
                
            end % END OF MW LOOP
            
            %% SPLIT TRIALS INTO IDX / NDX
            %  THEY ARE THEN SAVED IN VARIABLES LIKE "numRpplIdx". TRIAL IS
            %  A RANDOM EFFECT HERE. THIS PART IS SO THE VALUES FROM THE PREVIOUS STEP (WITHIN BUNDLE) ARE SAVED
            %  OVER BUNDLES
            
            if maxMW == 1 % DON'T TAKE ALL MW, INSTEAD GET THE MAXIMUM MW
                numRppl      = max(numRppl, [], 2);
                lenRppl      = max(lenRppl, [], 2);
                rpplDen      = max(rpplDen, [], 2);
                rpplSpkCoinc = max(rpplSpkCoinc, [], 2);
            end
           
            %% INDEXED
            %  INDEXED TRIALS - NUMBER OF RIPPLES
            tmp = numRppl( idxTrl,:);
            numRpplIdx = [numRpplIdx; tmp(:)]; % IF WE HAVE ALL 8 MW AT THIS POINT, THIS WILL POOL IT INTO ONE VECTOR
            
            %  INDEXED TRIALS - LENGTH OF RIPPLES
            tmp = lenRppl( idxTrl,:);
            lenRpplIdx = [lenRpplIdx; tmp(:)];
            
            %  INDEXED TRIALS - RIPPLE DENSITY
            tmp = rpplDen( idxTrl,:);
            rpplDenIdx = [rpplDenIdx; tmp(:)];
            
            %  INDEXED TRIALS - SPIKE RIPPLE COINCIDENCE
            tmp = rpplSpkCoinc( idxTrl,:);
            coincIdx = [coincIdx; tmp(:)];
            
            %% NON-INDEXED
            %  NON-INDEXED TRIALS - NUMBER OF RIPPLES
            tmp = numRppl(~idxTrl,:);
            numRpplNdx = [numRpplNdx; tmp(:)];
            
            %  NON-INDEXED TRIALS - LENGTH OF RIPPLES
            tmp = lenRppl(~idxTrl,:);
            lenRpplNdx = [lenRpplNdx; tmp(:)];
              
            %  NON-INDEXED TRIALS - RIPPLE DENSITY
            tmp = rpplDen(~idxTrl,:);
            rpplDenNdx = [rpplDenNdx; tmp(:)];
            
            %  NON-INDEXED TRIALS - SPIKE RIPPLE COINCIDENCE
            tmp = rpplSpkCoinc(~idxTrl,:);
            coincNdx = [coincNdx; tmp(:)];
        end
        
        %% VISUALIZE
        figure('units','normalized','outerposition',[0 0 1 1]);
        
        % FIGURE TITLE
        tits = '';
        if tw == 0
            tits = strjoin({tits, 'Whole Trial |'});
        elseif tw == 1
            tits = strjoin({tits, 'Pre-Cue |'});
        elseif tw == 2
            tits = strjoin({tits, 'Peri-Cue |'});
        elseif tw == 3
            tits = strjoin({tits, 'Peri-Resp |'});
        end
        
        if maxMW == 0
            tits = strjoin({tits, 'All Microwires'});
        else
            tits = strjoin({tits, 'Only Microwire with maximum value'});
        end
%         sgtitle(tits)
        
        % RIPPLE NUMBER
        subplot(211)
        [p, ~] = ranksum(numRpplIdx, numRpplNdx, 'tail', 'right');
        [p_perm] = perm_ranksum(numRpplIdx, numRpplNdx);
        title({sprintf('Number of Ripples (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
            sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(numRpplIdx), mean(numRpplIdx), median(numRpplIdx)), ...
            sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(numRpplNdx), mean(numRpplNdx), median(numRpplNdx))});
        hold on
        hand1 = histogram(numRpplIdx, 0:1:10, 'Normalization','probability');
        hand2 = histogram(numRpplNdx, 0:1:10, 'Normalization','probability');
        hand2.FaceColor = [1,0,0];
        plot([mean(numRpplIdx) mean(numRpplIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 4, 'linestyle', '-');
        plot([mean(numRpplNdx) mean(numRpplNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '-');
        plot([median(numRpplIdx) median(numRpplIdx)], get(gca, 'YLim'), 'color', 'b',     'linew', 4, 'linestyle', '--');
        plot([median(numRpplNdx) median(numRpplNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '--');
        legend('Indexed', 'Non-Indexed')
        xlabel(tits);
        
        % RIPPLE LENGTH
        subplot(212)
        [p, ~] = ranksum(lenRpplIdx, lenRpplNdx, 'tail', 'right');
        [p_perm] = perm_ranksum(lenRpplIdx, lenRpplNdx);
        title({sprintf('Length of Ripples (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
            sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(lenRpplIdx), nanmean(lenRpplIdx), nanmedian(lenRpplIdx)), ...
            sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(lenRpplNdx), nanmean(lenRpplNdx), nanmedian(lenRpplNdx))});
        
        hold on
        hand1 = histogram(lenRpplIdx, 30:5:120, 'Normalization','probability');
        hand2 = histogram(lenRpplNdx, 30:5:120, 'Normalization','probability');
        hand2.FaceColor = [1,0,0];
        plot([nanmean(lenRpplIdx) nanmean(lenRpplIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 4, 'linestyle', '-');
        plot([nanmean(lenRpplNdx) nanmean(lenRpplNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '-');
        plot([nanmedian(lenRpplIdx) nanmedian(lenRpplIdx)], get(gca, 'YLim'), 'color', 'b',     'linew', 3, 'linestyle', '--');
        plot([nanmedian(lenRpplNdx) nanmedian(lenRpplNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '--');
        legend('Indexed', 'Non-Indexed')
        
        % SAVE FIGURE
        saveTits = '';
        if tw == 0
            saveTits = strjoin({saveTits, 'WholeTrial_'});
        elseif tw == 1
            saveTits = strjoin({saveTits, 'PreCue_'});
        elseif tw == 2
            saveTits = strjoin({saveTits, 'PeriCue_'});
        elseif tw == 3
            saveTits = strjoin({saveTits, 'PeriResp_'});
        end
        
        if maxMW == 0
            saveTits = strjoin({saveTits, 'allMW'});
        elseif maxMW == 1
            saveTits = strjoin({saveTits, 'maxMW'});
        end
        
        saveas(gcf,['X:\Luca\visu\rpplNum_rpplLen_', saveTits], 'png' )
        
        % VISUALIZE RIPPLE DENSITY
        close all
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(20,1,[4:20])
        
        [p, ~] = ranksum(rpplDenIdx, rpplDenNdx, 'tail', 'right');
        [p_perm] = perm_ranksum(lenRpplIdx, lenRpplNdx);
        title({sprintf('Density of Ripples (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
            sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(rpplDenIdx), mean(rpplDenIdx), median(rpplDenIdx)), ...
            sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(rpplDenNdx), mean(rpplDenNdx), median(rpplDenNdx))});
        hold on
        
        hand1 = histogram(rpplDenIdx, 0:0.0025:0.1, 'Normalization','probability');
        hand2 = histogram(rpplDenNdx, 0:0.0025:0.1, 'Normalization','probability');
        hand2.FaceColor = [1,0,0];
        
        plot([mean(rpplDenIdx) mean(rpplDenIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 2, 'linestyle', '-');
        plot([mean(rpplDenIdx) mean(rpplDenIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 2, 'linestyle', '-');
        plot([mean(rpplDenNdx) mean(rpplDenNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 2, 'linestyle', '-');
        
        plot([median(rpplDenIdx) median(rpplDenIdx)], get(gca, 'YLim'), 'color', 'b',     'linew', 2, 'linestyle', '--');
        plot([median(rpplDenNdx) median(rpplDenNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 2, 'linestyle', '--');
        
        legend('Indexed', 'Non-Indexed')
        %         sgtitle(tits)
        xlabel(tits)
        
        saveas(gcf,['X:\Luca\visu\rpplDen_', saveTits], 'png' )
        close all
        
        % VISUALIZE RIPPLE COINCIDENCES
        close all
        figure('units','normalized','outerposition',[0 0 1 1]);
                
        [p, ~] = ranksum(coincIdx, coincNdx, 'tail', 'right');
        [p_perm] = perm_ranksum(coincIdx, coincNdx);
        title({sprintf('Spike-Ripple Coincidences (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
            sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(coincIdx), mean(coincIdx), median(coincIdx)), ...
            sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(coincNdx), mean(coincNdx), median(coincNdx))});
        hold on
        hand1 = histogram(coincIdx, 0:1:10, 'Normalization','probability');
        hand2 = histogram(coincNdx, 0:1:10, 'Normalization','probability');
        hand2.FaceColor = [1,0,0];
        plot([mean(coincIdx) mean(coincIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 4, 'linestyle', '-');
        plot([mean(coincNdx) mean(coincNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '-');
        plot([median(coincIdx) median(coincIdx)], get(gca, 'YLim'), 'color', 'b',     'linew', 4, 'linestyle', '--');
        plot([median(coincNdx) median(coincNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '--');
        legend('Indexed', 'Non-Indexed')
        %         sgtitle(tits)
        xlabel(tits);
        
        saveas(gcf,['X:\Luca\visu\rpplCoin_', saveTits], 'png' )
        close all

    end
end