function findIndexSU_varWind(plvl, fixTH)
for encWindow = 5 % 1:22
%     if encWindow~=10
%         continue
%     end
    for retWindow = 7 % 1:22
%         if retWindow~=7
%             continue
%         end
        
        clearvars -except plvl fixTH encWindow retWindow
        [timeWindowENC, timeWindowRET, lockingENC, lockingRET, onlyCues] = findIndexSU_timeWindows(encWindow, retWindow);
        
        allDatPath = 'X:\Luca\data\allSbj';
        allDat = dir([allDatPath, filesep, 'allSU*']);
        suNum = 0;
        encVec_allSU = [];
        retVec_allSU = [];
        numperm = 10000;
        empNumFP2 = 0;
        savingpath = ['X:\Luca\indexSUvisu\\varWindow_ENC_', num2str(timeWindowENC(1)), num2str(timeWindowENC(2)), num2str(lockingENC(1)), num2str(lockingENC(2)), '_RET_', num2str(timeWindowRET(1)), num2str(timeWindowRET(2)), num2str(lockingRET(1)), num2str(lockingRET(2)), '_DP_p',num2str(plvl),'_th',num2str(fixTH)];
%         if exist(savingpath, 'dir') == 7
%             continue
%         end
        
        mkdir(savingpath)
        
        for ia = 1:size(allDat,1)
            
            % load allSU file that contains the waveshape and spike timestamps
            load([allDat(ia).folder, filesep, allDat(ia).name])
            if isempty(allSU)
                break
            end
            
            % get bidsID  & sesh number
            bidsID = allDat(ia).name;
            bidsID([1:6 end-6:end]) = [];
            subj = sub_ID_conversion(bidsID, 'yes');
            sesh = allDat(ia).name;
            sesh = sesh(end-5:end-4);
            
            %                 subjID = [subjID, '_', sesh];
            
            if strcmp(sesh, '1b')
                subj= 'P07ERL';
                sesh = 'S1b';
            end
            
            cd(['X:/Luca/data', filesep, subj, filesep, sesh])
            abc = dir; cd(abc(3).name);
            p2d = [cd, filesep];
            
            %% load logfile
            [tableTemplate, hitsIdx, ~, ~, ~, retTrigger, encTrigger, ~, ~] = loadLogs(p2d, 1);
            
            spiketimes_segmentedENC = {};
            spiketimes_segmentedRET = {};
            for su = 1:size(allSU,1)
                [spiketimes_segmentedENC(su,:), twlengthEnc] = insertSpiketimes2(encTrigger, allSU{su,3}./32000, lockingENC, timeWindowENC);
                [spiketimes_segmentedRET(su,:), twlengthRet] = insertSpiketimes2(retTrigger, allSU{su,3}./32000, lockingRET, timeWindowRET);
            end
            
            % transform from spiketimes to spikenumber
            spikeNumber_segmentedENC = cellfun(@length, spiketimes_segmentedENC) / abs(twlengthEnc);
            spikeNumber_segmentedRET = cellfun(@length, spiketimes_segmentedRET) / abs(twlengthRet);
            
            % extract hits / delete misses
            spikeNumber_segmentedENC = spikeNumber_segmentedENC (:,hitsIdx);
            spikeNumber_segmentedRET = spikeNumber_segmentedRET (:,hitsIdx);
            
            [normENC, normRET]= normSpikeNumber_double(spikeNumber_segmentedENC, spikeNumber_segmentedRET); % normalizes spikenumber (ENC + RET separate)
            numSU = size(normENC,1);
            
            %         loop over all single units in that session
            for su = 1 : size(normENC,1)
                suNum = suNum+1;
                wireName = allSU{su,1};
                % dot product of both TS before shuffling
                origDP  = normENC(su,:) .* normRET(su,:);
                
                origENC = normENC(su,:); % original encoding  trial series (presentation order)
                origRET = normRET(su,:); % original retrieval trial series (presentation order)
                
                %             logical for each stimulus category (ff/pp/fp)
                %             used for color coding the trials in the TS
                ffIdx{suNum}   = ~cellfun(@isempty, strfind(tableTemplate{1,:}, 'ff'));
                ppIdx{suNum}   = ~cellfun(@isempty, strfind(tableTemplate{1,:}, 'pp'));
                fpIdx{suNum}   = ~cellfun(@isempty, strfind(tableTemplate{1,:}, 'fp'));
                himiIdx        =  cellfun(@isempty, strfind(tableTemplate{2,:}, 'hit')); % if hitMiss is 'hit' himiIdx will have a 1 if a trial was a miss
                
                %             kick out the trials with either hits or miss (depending on hitMiss)
                ffIdx{suNum}(himiIdx) = [];
                ppIdx{suNum}(himiIdx) = [];
                fpIdx{suNum}(himiIdx) = [];
                
                permdxR = []; % delete the varible because different subjects have different trial lengths
                permdxE = [];
                for prm = 1:numperm
                    %                 for DP shuffle across categories
                    permdxR = randperm(size(hitsIdx,1)); % reshuffle retrieval
                    permdxE = randperm(size(hitsIdx,1)); % reshuffle encoding
                    permDP{suNum}(:,prm)  =  normENC(su,permdxE) .* normRET(su,permdxR); % permutated dot product between the two time series, note that here we have to shuffle encoding as well!
                    
                    %             %% enc + ret separate
                    permEnc{suNum}(:,prm) = normENC(su, permdxE);
                    permRet{suNum}(:,prm) = normRET(su, permdxR);
                end
                
                % calculate the threshold
                % DP
                dynTH                =  prctile(permDP{suNum}(:),plvl);
                threshDP(suNum)      = dynTH;
                sensTrls_perm{suNum} = permDP{suNum} >= threshDP(suNum) & abs(permEnc{suNum}) >= fixTH & abs(permRet{suNum}) >= fixTH;
                
                %%
                sensTrls = origDP >= threshDP(suNum) & origENC >= fixTH & origRET >= fixTH; % sensitive trials
                ffTrls = and(ffIdx{suNum}, sensTrls);
                ppTrls = and(ppIdx{suNum}, sensTrls);
                fpTrls = and(fpIdx{suNum}, sensTrls);
                
                if onlyCues == 0
                    putIdx = sum(ffTrls) >= 1 && sum(ppTrls) >=1;
                elseif onlyCues == 1
                    putIdx = sum([ffTrls, ppTrls, fpTrls]) >= 2;
                end
                % tests if a SU is a FP2 neuron
                if putIdx
                    %             allDP       = [allDP, origDP]; % note down their DP in the main diagnonal
                    empNumFP2   = empNumFP2+1; % count how many FP2 I found empirically
                    %             whichSUfp2  = [whichSUfp2, suNum]; % mark their number
                    
                    %% % extract the timestamps between cue onset and 1s post cue firing
                    idxTrls = logical(sum([ffTrls; ppTrls; fpTrls]));
                    encVec = vertcat(spiketimes_segmentedENC{su, idxTrls});
                    encVec(or(encVec<0, encVec>1)) = [];
                    
                    retVec = vertcat(spiketimes_segmentedRET{su, idxTrls});
                    retVec(or(retVec<0, retVec>1)) = [];
                    
                    if and(~isempty(encVec), ~isempty(retVec)) %
                        encVec_allSU = [encVec_allSU, {encVec}];
                        retVec_allSU = [retVec_allSU, {retVec}];
                    end
                    
                    
                    %%                     visualisation
                    mhandle = figure;
                    subplot(2,1,1)
                    hold on
                    
                    %                     plotting
                    plot(normENC(su,:), 'linew',3)
                    plot(normRET(su,:), 'linew',3)
                    
                    %                     aesthetics
                    xlim([1 length(ffIdx{suNum})])
                    %             xlabel('Trial (ordered by encoding presentation');
                    hand = legend({'ENC', 'RET', ''});
                    ylabel('Standardized firing rate');
                    title([bidsID, '-', sesh, sprintf('-SU%.0f', su)]);
                    set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
                    axH = gca;
                    axH.YAxis.FontWeight = 'bold';
                    axH.XAxis.FontWeight = 'bold';
                    axH.FontSize = 24;
                    
                    %% DP
                    subplot(2,1,2)
                    hold on
                    plot(origDP, 'linew',3)
                    plot([1 length(ffIdx{suNum})], [threshDP(suNum) threshDP(suNum)], 'b--', 'linew', 2)
                    
                    recPos = get(gca,'Ylim');
                    startP = -0.5;
                    
                    % this marks only the trials in which the dp exceeds threshold
                    for iy = 1 : length(ffIdx{suNum})
                        startP = startP+1;
                        if ffIdx{suNum}(iy)==1 && sensTrls(iy) ==1
                            r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[0 0 0 0.4]); % dark
                        elseif ppIdx{suNum}(iy)==1 && sensTrls(iy) ==1
                            r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[1 1 1 0.35]); % white
                        end
                    end
                    hold off
                    
                    % aesthetics
                    xlim([1 length(ffIdx{suNum})])
                    xlabel('Trial (ordered by encoding presentation)');
                    ylabel('DP ');
                    set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
                    axH = gca;
                    axH.YAxis.FontWeight = 'bold';
                    axH.XAxis.FontWeight = 'bold';
                    axH.FontSize = 24;
                    %                     saving figure
                    cd(savingpath)
                    saveas(gcf, [bidsID, '-', sesh, sprintf('-SU%.0f',su)], 'png');
                    
                    %
                    sensTrlsFFPP = or(ffIdx{suNum}==1 & sensTrls ==1, ppIdx{suNum}==1 & sensTrls ==1);
                    save( sprintf('SU_%d.mat',size(dir(('SU_*.mat')),1)+1), 'bidsID', 'sesh', 'wireName', 'su', 'sensTrlsFFPP', 'ffTrls', 'ppTrls')
                    
                end
                
            end
        end
        
        %% see how many FP2 I can expect under the null
        % ffIdx{numSU}         - a cell with a 1xTRL logical for each SU. The logical indexes which trials were ff trials
        % sensTrls_perm{numSU} - a cell with a TRLx10.000 logical for each SU. The logical indexes to which trial the SU is sensitive to in each of the 10.000 permutations
        
        %     load('X:\Luca\anCueNames.mat', 'anCueNames');
        %     nperm = 100000;
        %
        %     simil = zeros(1,nperm);
        %     for ia=1:nperm
        %         idx = randperm(size(anCueNames,1),2);
        %         [vec1, vec2] = anCueNames{idx,3};
        %         simil(ia) = cosSimil(vec1,vec2);
        %     end
        % killPermIU = sum(simil>conceptPermTH)/nperm; % what are the chances two images randomly chosen from our distribution of animals are at least as similar as the least similar concept neuron that we kicked out
        
        permNumFP2_glob = zeros([1 numperm]);
        %         x = zeros(1,10000);
        %         x(1:killPermIU*10000) = 1;
        for ih = 1:numperm % do this 10.000 times
            for numSU = 1:size(sensTrls_perm,2)
                perm = randi(numperm); % this actually might not be necessary
                
                ff = sum(ffIdx{numSU}(sensTrls_perm{numSU}(:, perm))); % take out a random permutation from that SU
                pp = sum(ppIdx{numSU}(sensTrls_perm{numSU}(:, perm)));
                fp = sum(fpIdx{numSU}(sensTrls_perm{numSU}(:, perm)));
                %          && dpCor (rownum,perm) > threshDPcor(rownum) % threshDPcor would
                %          be a lower threshold because its not WC
                
                if onlyCues == 0
                    putIdx = sum(ff) >= 1 && sum(pp) >=1;
                elseif onlyCues == 1
                    putIdx = sum([ff, pp, fp]) >= 2;
                end
                % tests if a SU is a FP2 neuron
                if putIdx
                    %                 if ff > 0 && pp > 0
                    %             %% chance that it is a concept neuron
                    %             x = x(randperm(length(x)));
                    %             if x(1) == 1
                    %                 continue
                    %             end
                    %%
                    permNumFP2_glob(ih) = permNumFP2_glob(ih) + 1;
                end
            end
        end
        th_fp2_glob = prctile(permNumFP2_glob,95) % gives significance / th at 15
        sumnum = sum(empNumFP2>=permNumFP2_glob) % gives 9996
        empNumFP2 % at 21
        cd(savingpath)
        %         save('globalTest', 'th_fp2_glob', 'sumnum', 'empNumFP2', 'encVec_allSU', 'retVec_allSU', 'permNumFP2_glob')
        save('globalTest', 'th_fp2_glob', 'sumnum', 'empNumFP2','encVec_allSU', 'retVec_allSU', 'permNumFP2_glob', 'ffIdx', 'ppIdx', 'fpIdx', 'sensTrls_perm')
        
        close all
    end
end
end % end of function