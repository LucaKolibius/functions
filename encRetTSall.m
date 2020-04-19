%% permutation testing
% all SU time series between encoding and retrieval
clear
close all
fetch_allSubj
suNum = 0;
numperm = 10000;
permCor = zeros(377, numperm); % preallocation for the permutated correlation. x single units, y permutations
% permDP = zeros(377, numperm);  % preallocation for the permutated dot product. x single units, y permutations
empNumSigSU = 0;
empNumFP2  = 0;
plvl(1) = 95; plvl(2) = 95; % plvl für correlation/dp
hitMiss = 'hit'; % hit / miss
fixTH = 2;
allDP = [];
allxDiagDP = [];
noFP2_DP = [];
noFP2_xDiagDP = [];
noCor_DP = [];
noCor_xDiagDP = [];
whichSUfp2 = [];
permNumFP2 = [];
permDP = {};

for sbj = 1 : size(allSubj,1)
    subjID = allSubj{sbj};
    disp(subjID);
    try
        cd X:/Luca_old/data
    catch
        cd /media/ldk898/rds-share/Luca/data
    end
    
    mSubject = subjID(1:end-3);
    mSession = subjID(end-1:end);
    
    % if the session name is called 1b then this line prevents an error during cd
    mSubject(regexp(mSubject,'_')) = [];
    if isempty(regexp(mSession,'S', 'ONCE'))
        mSession = ['S', mSession];
    end
    
    cd(mSubject)
    cd(mSession)
    abc = dir;
    cd(abc(3).name)
    cd advancedAnalysis\elecLoc\spikenumbers
    
    % load tables for enc and ret
    abc = dir('spikeNumber_hipp_*');
    load(abc.name, 'encSpikeNumber_cueLocked_h', 'retSpikeNumber_cueLocked_h')
    
    allTrials = size(encSpikeNumber_cueLocked_h,2);
    if size(encSpikeNumber_cueLocked_h,1)>2
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell2(allTrials, encSpikeNumber_cueLocked_h, retSpikeNumber_cueLocked_h, hitMiss); % extracts the absolute spikenumbers for encoding and retrieval. CAREFUL: only hits and ordered by presentation time at encoding
    else
        continue
    end
    
    if size(tempCell_enc_h,2)>1 % in case there are no misses
        [~ , normENC, normRET] = normSpikeNumber(tempCell_enc_h, tempCell_ret_h); % normalizes spikenumber and correlates enc+ret for RSA
    else
        continue
    end
    subjID(regexp(subjID,'_')) = '-';
    numSU = size(normENC,1);
    
    % loop over all single units in that session
    for su = 1 : size(normENC,1)
        suNum = suNum+1;
        
        % correlate normalized encoding trial series with normalized retrieval trial series (before shuffling)
        tempCor = corrcoef(normENC(su,:), normRET(su,:));
        origCor = tempCor(2); 
        
        % dot product of both TS before shuffling
        origDP  = normENC(su,:) .* normRET(su,:);
        
        % dot product of the cross diagonal of the RSA
        xDiagDP = mXDiagDP(encSpikeNumber_cueLocked_h, retSpikeNumber_cueLocked_h, hitMiss, su);
        
        origENC = normENC(su,:); % original encoding  trial series (presentation order)
        origRET = normRET(su,:); % original retrieval trial series (presentation order)
                
        % logical for each stimulus category (ff/pp/fp)
        % used for color coding the trials in the TS
        ffIdx{suNum}   = ~cellfun(@isempty, strfind(encSpikeNumber_cueLocked_h{1,:}, 'ff'));
        ppIdx{suNum}   = ~cellfun(@isempty, strfind(encSpikeNumber_cueLocked_h{1,:}, 'pp'));
        fpIdx{suNum}   = ~cellfun(@isempty, strfind(encSpikeNumber_cueLocked_h{1,:}, 'fp'));
        himiIdx        =  cellfun(@isempty, strfind(encSpikeNumber_cueLocked_h{2,:}, hitMiss)); % if hitMiss is 'hit' himiIdx will have a 1 if a trial was a miss
        
        % kick out the trials with either hits or miss (depending on hitMiss)
        ffIdx{suNum}(himiIdx) = [];
        ppIdx{suNum}(himiIdx) = [];
        fpIdx{suNum}(himiIdx) = [];
        
        % ffIdx2 is an index of which trials are hits and ff
        % this is for the WC permutation
        ffIdx2{suNum}  =  find(ffIdx{suNum}==1);
        ppIdx2{suNum}  =  find(ppIdx{suNum}==1);
        fpIdx2{suNum}  =  find(fpIdx{suNum}==1);
        
        permdxR = []; % delete the varible because different subjects have different trial lengths
        permdxE = [];
        for prm = 1:numperm
            % permutate within one category. ffIdx{rownum} is the index of when an
            % ff trial occurs for accapted SU 'rownum'. This could be e.g.
            % 4 31 38 2 9 13 22 36 40 16 19 37 28
            % Then we randomly permutate this index WITHIN each category.
            % This results in ffRP etc.
            ffRP            =  ffIdx2{suNum}(randperm(size(ffIdx2{suNum},2), size(ffIdx2{suNum},2)));
            ppRP            =  ppIdx2{suNum}(randperm(size(ppIdx2{suNum},2), size(ppIdx2{suNum},2)));
            fpRP            =  fpIdx2{suNum}(randperm(size(fpIdx2{suNum},2), size(fpIdx2{suNum},2)));
            
            % Because I do not want to simply merge these together, I then
            % insert this randomized WC index back into the old ff/pp/fp
            % position
            permdxR(ffIdx{suNum}) = ffRP;
            permdxR(ppIdx{suNum}) = ppRP;
            permdxR(fpIdx{suNum}) = fpRP;
            
            permdxE(ffIdx{suNum}) = ffRP;
            permdxE(ppIdx{suNum}) = ppRP;
            permdxE(fpIdx{suNum}) = fpRP;
            
             % recalculate correlation with shuffled label for retrieval
            tempCor               =  corrcoef(normENC(su,:),  normRET(su,permdxR));
            permCor(suNum,prm)    =  tempCor(2); % extract correlation
            
            % for DP shuffle across categories
            permdxR = randperm(size(tempCell_ret_h,2)); % reshuffle retrieval
            permdxE = randperm(size(tempCell_ret_h,2)); % reshuffle encoding
            permDP{suNum}(:,prm)  =  normENC(su,permdxE) .* normRET(su,permdxR); % permutated dot product between the two time series, note that here we have to shuffle encoding as well!
            
            tmpCor = corrcoef(normENC(su,:), normRET(su,permdxR));
            dpCor(suNum, prm) = tmpCor(2);
            %             permENC{rownum}(:,prm) = normENC(su,permdxE); % save each iteration of shuffled encoding trial series for SU rownum
            %             permRET{rownum}(:,prm) = normRET(su,permdxR); % same for retrieval
        end
        
        % get rid of 1/-1 rows and NaNs
        if sum(isnan(permCor(suNum,:))) == numperm || (sum(permCor(suNum,:)==1 | permCor(suNum,:)==-1) == numperm)% if all permutations of that SU are nans or either -1 or 1
            % delete entries from that SU
            permCor(suNum,:) = [];
            dpCor(suNum,:) = [];
            permDP = permDP(1:end-1);
            %             permENC = permENC(1:end-1);
            %             permRET = permRET(1:end-1);
            ffIdx = ffIdx(1:end-1);
            ppIdx = ppIdx(1:end-1);
            fpIdx = fpIdx(1:end-1);
            
            % decrease SU counter by 1
            suNum = suNum-1;
            continue
        else % if I keep the permutation for that SU, calculate the threshold
            thresh       (suNum)  =  prctile( permCor(suNum,:),  plvl(1)); % somehow does not give the exact same result as above.. at least not for sbj= 2, su=3;
            dynTH                 =  prctile(permDP{suNum}(:),plvl(2));
            threshDP     (suNum)  =  max(dynTH, fixTH);
            
            
%             threshDPcor  (rownum)  =  prctile(dpCor (rownum,:), plvl(1));
            %             threshENC (rownum)  =  prctile( permENC{rownum}(:), plvl(2));
            %             threshRET (rownum)  =  prctile( permRET{rownum}(:), plvl(2));
            %             permCorSort      = sort(permCor(rownum,:));
            %             thresh(rownum)   = permCorSort(numperm-numperm*plvl);
        end
        
        %% find out how often a single! SU is FP2 (FP squared / FFPP)
                % at this stage I do not take FP2 into account. Just look if the shuffled SU
                % is sensitive toward the same trial at encoding and retrieval
                %         sensTrl_permEnc = permENC{rownum} >= threshENC(rownum);          % permutated SU sensitive to which encoding trials?
                %         sensTrl_permRet = permRET{rownum} >= threshRET(rownum);          % same for retrieval
                %         sensTrls_perm{rownum}   = sensTrl_permEnc+sensTrl_permRet == 2;  % SU sensitive to same encoding and retrieval trial?
                sensTrls_perm{suNum} = permDP{suNum} >= threshDP(suNum);
                
                % counter how often I found that this shuffled SU is sensitive towards FP2 at encoding and retrieval
                counter = 0;
                for col = 1:numperm
                    ffTrls_perm = sum(ffIdx{suNum}(sensTrls_perm{suNum}(:,col)));
                    ppTrls_perm = sum(ppIdx{suNum}(sensTrls_perm{suNum}(:,col)));
                    
                    if ffTrls_perm >0 && ppTrls_perm >0
                        counter = counter+1; % out of the 10.000 permutations, how often is this SU sensitive to FP2?
                    end
                end
                permNumFP2 = [permNumFP2, counter];
        
        %% if correlation is significant
        if origCor >= thresh(suNum) % if the correlation between encoding and retrieval in the original trial series is significant
            empNumSigSU = empNumSigSU+1;
            whichSUfp2  = [whichSUfp2, suNum];
            
            
            %             sensTrl_enc = origENC >= threshENC;         % original SU is sensitive towards these encoding  trials
            %             sensTrl_ret = origRET >= threshRET;         % original SU is sensitive towards these retrieval trials
            %             sensTrls    = sensTrl_enc + sensTrl_ret == 2;  % SU is sensitive towards this trial at ENC + at RET
            
            sensTrls = origDP >= threshDP(suNum); % sensitive trials
            ffTrls = sum(ffIdx{suNum}(sensTrls));
            ppTrls = sum(ppIdx{suNum}(sensTrls));
            
            % tests if a SU is a FP2 neuron
            if ffTrls>1 && ppTrls >1 
                allDP       = [allDP, origDP]; % note down their DP in the main diagnonal
                allxDiagDP  = [allxDiagDP, xDiagDP]; % note down their DP in the cross diagnonal
                
                empNumFP2 = empNumFP2+1; % count how many FP2 I found empirically
                whichSUfp2 = [whichSUfp2, suNum]; % mark their number


                
                %                 % visualisation
                %                 mhandle = figure;
                %                 hold on
                %
                %                 % plotting
                %                 plot(normENC(su,:), 'linew',3)
                %                 plot(normRET(su,:), 'linew',3)
                %
                %                 % aesthetics
                %                 xlim([1 length(ffIdx)])
                %                 xlabel('Trial (ordered by encoding presentation');
                %                 hand = legend({'ENC', 'RET', ''});
                %                 ylabel('Standardized firing rate');
                %                 title([subjID, sprintf('-SU%.0f-corr%.2f', su, origCor)]);
                %                 set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
                %                 axH = gca;
                %                 axH.YAxis.FontWeight = 'bold';
                %                 axH.XAxis.FontWeight = 'bold';
                %                 axH.FontSize = 24;
                %
                %                 recPos = get(gca,'Ylim');
                %                 startP = -0.5;
                
                % this marks all trials
                %                 for iy = 1 : length(ffIdx)
                %                     startP = startP+1;
                %                     if ffIdx(iy)==1
                %                         r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[0 0 0 0.4]);
                %                     elseif ppIdx(iy)==1
                %                         r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[1 1 1 0.35]);
                %                     elseif fpIdx(iy)==1
                %                         r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[0 0 0 0.2]); % [startP -1 1 6]
                %                     end
                %                 end
                %                 hold off
                
                % % this marks only the trials in which the dp exceeds threshold
                % for iy = 1 : length(ffIdx)
                %     startP = startP+1;
                %     if ffIdx(iy)==1 && sensTrls(iy) ==1
                %         r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[0 0 0 0.4]); % dark
                %     elseif ppIdx(iy)==1 && sensTrls(iy) ==1
                %         r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[1 1 1 0.35]); % white
                %     elseif fpIdx(iy)==1 && sensTrls(iy) ==1
                %         r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[0 0 0 0.2]); % [startP -1 1 6] % light grey
                %     end
                % end
                % hold off
                
                %                 % saving figure
                %                 cd X:\Luca\data\SUvisu\indexSU
                %                 saveas(gcf, [subjID, sprintf('-SU%.0f',su)], 'png');
            
            else % if the correlation exceeds threshold, but it is not a FP2 unit
                noFP2_DP = [noFP2_DP, origDP];
                noFP2_xDiagDP = [noFP2_xDiagDP, xDiagDP];
            end
        else % if the correlation does not exceed the threshold (FP2 is not even tested then)
            noCor_DP       = [noCor_DP, origDP];
            noCor_xDiagDP  = [noCor_xDiagDP, xDiagDP];
        end
    end
end

%% How many significantly correlating SU do I expect under the H0?
% extract randomly one permutation from each SU
for ih = 1:numperm
    permdxASU = randperm(10000,size(permCor,1))';
    temp = permCor([1:size(permCor,1)], permdxASU);
    permvec(:,ih) = temp(logical(eye(size(permCor,1)))); % a vector of randomly drawn values from all SU (one value per SU)
    numsign(ih) = sum(permvec(:,ih)>=thresh'); % the number of randomly drawn h0 correlations that are signficiantly above threshold
end
th_cor_glob = prctile(numsign,95);
sum(empNumSigSU>numsign)
max(numsign)

% numsignSort = sort(numsign);
% numsignSort(numperm-numperm*0.05) % threshold
% numsignSort(end)


%% see how many FP2 I can expect under the null
% ffIdx{numSU}         - a cell with a 1xTRL logical for each SU. The logical indexes which trials were ff trials
% sensTrls_perm{numSU} - a cell with a TRLx10.000 logical for each SU. The logical indexes to which trial the SU is sensitive to in each of the 10.000 permutations

permNumFP2_glob = zeros([1 numperm]);
for ih = 1:numperm % do this 10.000 times
    for numSU = 1:size(sensTrls_perm,2)
        perm = randi(numperm); % this actually might not be necessary
        
        ff = sum(ffIdx{numSU}(sensTrls_perm{numSU}(:, perm))); % take out a random permutation from that SU
        pp = sum(ppIdx{numSU}(sensTrls_perm{numSU}(:, perm)));
        
%          && dpCor (rownum,perm) > threshDPcor(rownum) % threshDPcor would
%          be a lower threshold because its not WC
        if ff > 0 && pp > 0 && dpCor (numSU,perm) > thresh(numSU) % && check if it correlates significantly between encoding and retrieval
            permNumFP2_glob(ih) = permNumFP2_glob(ih) + 1;
        end
    end
end
th_fp2_glob = prctile(permNumFP2_glob,95) % gives significance / th at 15
sum(empNumFP2>=permNumFP2_glob) % gives 9996
empNumFP2 % at 21