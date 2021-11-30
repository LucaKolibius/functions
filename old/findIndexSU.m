% permutation testing
% all SU time series between encoding and retrieval

% for DP encRet = 0
% for encoding and retrieval individually encRet = 1;
function findIndexSU(plvl, fixTH, conceptPermTH)
close all
fetch_allSubj
suNum = 0;
numperm = 10000;
% permDP = zeros(377, numperm);  % preallocation for the permutated dot product. x single units, y permutations
empNumFP2  = 0;
% plvl = 95; % plvl für dp
hitMiss = 'hit'; % hit / miss
% fixTH = 3;
allDP = [];
noFP2_DP = [];
whichSUfp2 = [];
permNumFP2 = [];
permDP = {};
encVec_allSU = [];
retVec_allSU = [];

% encRet = 0; % 1 for enc+ret, 0 for DP
% switch encRet
%     case 0 % for DP
        mkdir([ 'X:\Luca\indexSUvisu\DP_p',num2str(plvl),'_th',num2str(fixTH) ])
        savingpath = ['X:\Luca\indexSUvisu\DP_p',num2str(plvl),'_th',num2str(fixTH)];
%     case 1 % for enc+ret
%         mkdir([ 'X:\Luca\indexSUvisu\encRet_p',num2str(plvl),'_th',num2str(fixTH) ])
%         savingpath = ['X:\Luca\indexSUvisu\encRet_p',num2str(plvl),'_th',num2str(fixTH)];
% end

for sbj = 1 : size(allSubj,1)
    subjID = allSubj{sbj};
    disp(subjID);
    try
        cd X:/Luca/data
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
    
    cd ..\spiketimes
    abc = dir('spiketimes_hipp_*');
    load(abc.name, 'encSpiketimes_cueLocked_h', 'retSpiketimes_cueLocked_h');
    
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
    
    %         loop over all single units in that session
    for su = 1 : size(normENC,1)
        suNum = suNum+1;
                
        % dot product of both TS before shuffling
        origDP  = normENC(su,:) .* normRET(su,:);
                
        origENC = normENC(su,:); % original encoding  trial series (presentation order)
        origRET = normRET(su,:); % original retrieval trial series (presentation order)
        
        %             logical for each stimulus category (ff/pp/fp)
        %             used for color coding the trials in the TS
        ffIdx{suNum}   = ~cellfun(@isempty, strfind(encSpikeNumber_cueLocked_h{1,:}, 'ff'));
        ppIdx{suNum}   = ~cellfun(@isempty, strfind(encSpikeNumber_cueLocked_h{1,:}, 'pp'));
        fpIdx{suNum}   = ~cellfun(@isempty, strfind(encSpikeNumber_cueLocked_h{1,:}, 'fp'));
        himiIdx        =  cellfun(@isempty, strfind(encSpikeNumber_cueLocked_h{2,:}, hitMiss)); % if hitMiss is 'hit' himiIdx will have a 1 if a trial was a miss
        
        %             kick out the trials with either hits or miss (depending on hitMiss)
        ffIdx{suNum}(himiIdx) = [];
        ppIdx{suNum}(himiIdx) = [];
        fpIdx{suNum}(himiIdx) = [];
        
        %             ffIdx2 is an index of which trials are hits and ff
        %             this is for the WC permutation
        ffIdx2{suNum}  =  find(ffIdx{suNum}==1);
        ppIdx2{suNum}  =  find(ppIdx{suNum}==1);
        fpIdx2{suNum}  =  find(fpIdx{suNum}==1);
        
        permdxR = []; % delete the varible because different subjects have different trial lengths
        permdxE = [];
        
        for prm = 1:numperm
%             %                 permutate within one category. ffIdx{rownum} is the index of when an
%             %                 ff trial occurs for accapted SU 'rownum'. This could be e.g.
%             %                 4 31 38 2 9 13 22 36 40 16 19 37 28
%             %                 Then we randomly permutate this index WITHIN each category.
%             %                 This results in ffRP etc.
%             ffRP            =  ffIdx2{suNum}( randperm( size(ffIdx2{suNum},2), size(ffIdx2{suNum},2) ));
%             ppRP            =  ppIdx2{suNum}( randperm( size(ppIdx2{suNum},2), size(ppIdx2{suNum},2) ));
%             fpRP            =  fpIdx2{suNum}( randperm( size(fpIdx2{suNum},2), size(fpIdx2{suNum},2) ));
%             
%             %                 Because I do not want to simply merge these together, I then
%             %                 insert this randomized WC index back into the old ff/pp/fp
%             %                 position
%             permdxR(ffIdx{suNum}) = ffRP;
%             permdxR(ppIdx{suNum}) = ppRP;
%             permdxR(fpIdx{suNum}) = fpRP;
            
            %                 for DP shuffle across categories
            permdxR = randperm(size(tempCell_ret_h,2)); % reshuffle retrieval
            permdxE = randperm(size(tempCell_ret_h,2)); % reshuffle encoding
            permDP{suNum}(:,prm)  =  normENC(su,permdxE) .* normRET(su,permdxR); % permutated dot product between the two time series, note that here we have to shuffle encoding as well!
        
%             %% enc + ret separate
% permuting this doesn't make any sense because its just the trial series
% vector ( max([1 2 3 4 5]) == max([3 5 1 2 4]) )
            permEnc{suNum}(:,prm) = normENC(su, permdxE);
            permRet{suNum}(:,prm) = normRET(su, permdxR);        
        end
        
        % calculate the threshold
        % DP
        dynTH             =  prctile(permDP{suNum}(:),plvl);
%         threshDP (suNum)  =  max(dynTH, fixTH);
threshDP(suNum) = dynTH;
        
%         % ENC + RET
%         dynTH_enc = prctile(permEnc{suNum}(:), plvl);
%         dynTH_ret = prctile(permRet{suNum}(:), plvl);
%         
%         threshEnc {suNum} = max(dynTH_enc, fixTH);
%         threshRet {suNum} = max(dynTH_ret, fixTH);
        
        %                 % find out how often a single! SU is FP2 (FP squared / FFPP)
        %                 at this stage I do not take FP2 into account. Just look if the shuffled SU
        %                 is sensitive toward the same trial at encoding and retrieval
        
%         if encRet == 0
            sensTrls_perm{suNum} = permDP{suNum} >= threshDP(suNum) & abs(permEnc{suNum}) >= fixTH & abs(permRet{suNum}) >= fixTH;
%         elseif encRet == 1
%             sensTrls_perm{suNum} = and(permEnc{suNum} >= threshEnc{suNum}, permRet{suNum}>=threshRet{suNum});
%         end
        
%         % %                 counter how often I found that this shuffled SU is sensitive towards FP2 at encoding and retrieval
%         counter = 0;
%         for col = 1:numperm
%             % the permuted data obv dont conform to the ffIdx and ppIdx
%             % labels anymore. But that shouldnt matter here, because of the
%             % scrambling. The important bit is that the number of ff and pp
%             % trls is taken into account
%             ffTrls_perm = sum(ffIdx{suNum}(sensTrls_perm{suNum}(:,col)));
%             ppTrls_perm = sum(ppIdx{suNum}(sensTrls_perm{suNum}(:,col)));
%             
%             if ffTrls_perm >0 && ppTrls_perm >0
%                 counter = counter+1; % out of the 10.000 permutations, how often is this SU sensitive to FP2?
%             end
%         end
%         permNumFP2 = [permNumFP2, counter];
        
       
        %%
        sensTrls = origDP >= threshDP(suNum) & origENC >= fixTH & origRET >= fixTH; % sensitive trials
        ffTrls = and(ffIdx{suNum}, sensTrls);
        ppTrls = and(ppIdx{suNum}, sensTrls);
        
%         if encRet == 1
%             sensTrls_enc = origENC >= threshEnc{suNum}; % sensitive trials for encoding  series
%             sensTrls_ret = origRET >= threshRet{suNum}; % sensitive trials for retrieval series
%             
%             ffTrls = and( and(ffIdx{suNum},sensTrls_enc), and(ffIdx{suNum},sensTrls_ret) ); % gives you the FF trials that the SU is sensitive towards at encoding AND retrieval
%             ppTrls = and( and(ppIdx{suNum},sensTrls_enc), and(ppIdx{suNum},sensTrls_ret) ); % same for PP trials
%         end
        
        % tests if a SU is a FP2 neuron
        if sum(ffTrls) >= 1 && sum(ppTrls) >=1
            allDP       = [allDP, origDP]; % note down their DP in the main diagnonal
            empNumFP2   = empNumFP2+1; % count how many FP2 I found empirically
            whichSUfp2  = [whichSUfp2, suNum]; % mark their number          
            
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
            title([subjID, sprintf('-SU%.0f', su)]);
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
            xlabel('Trial (ordered by encoding presentation');
            ylabel('DP ');
            set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
            axH = gca;
            axH.YAxis.FontWeight = 'bold';
            axH.XAxis.FontWeight = 'bold';
            axH.FontSize = 24;
            %                     saving figure
            cd(savingpath)
            saveas(gcf, [subjID, sprintf('-SU%.0f',su)], 'png');
            
            %
            sensTrlsFFPP = or(ffIdx{suNum}==1 & sensTrls ==1, ppIdx{suNum}==1 & sensTrls ==1);
            save( sprintf('SU_%d.mat',size(dir(('SU*')),1)+1), 'subjID', 'su', 'sensTrlsFFPP', 'ffTrls', 'ppTrls')
            
            %% new vector with timestamps
            indexTrials = ffTrls+ppTrls; % which trials are indexed by this SU
            indexTrials = find(indexTrials == 1);
            for ix = 1:length(indexTrials)
                % encoding
                encVec = cell2mat(encSpiketimes_cueLocked_h{su+2,indexTrials(ix)}); % extract spiketimes
                encVec(or(encVec<0, encVec>1)) = []; % only consider spiketimes from cue onset to 1s post cue onset
                
                % repeat for retrieval
                retVec = cell2mat(retSpiketimes_cueLocked_h{su+2,indexTrials(ix)});
                retVec(or(retVec<0, retVec>1)) = [];
                
                % save in bigger varible (for all trials and SU)
                % but only if the trial has spikes in encoding AND
                % retrieval between cue onset and cue+1s
                if and(~isempty(encVec), ~isempty(retVec)) % 
                    encVec_allSU = [encVec_allSU, {encVec}]; 
                    retVec_allSU = [retVec_allSU, {retVec}];
                end
            end
%             encVec = [cell2mat(encSpiketimes_cueLocked_h{su+2, ffTrls}')]; % +2 because the first line is stimulus category and the second is hit/miss
%             retVec = [cell2mat(retSpiketimes_cueLocked_h{su+2, ffTrls}')];
      
            %% rasterplot
            
        else % if the correlation exceeds threshold, but it is not a FP2 unit
            noFP2_DP = [noFP2_DP, origDP];
        end
        
    end
end


%% see how many FP2 I can expect under the null
% ffIdx{numSU}         - a cell with a 1xTRL logical for each SU. The logical indexes which trials were ff trials
% sensTrls_perm{numSU} - a cell with a TRLx10.000 logical for each SU. The logical indexes to which trial the SU is sensitive to in each of the 10.000 permutations
    load('X:\Luca\anCueNames.mat', 'anCueNames');
    nperm = 100000;
    
    simil = zeros(1,nperm);
    for ia=1:nperm
        idx = randperm(size(anCueNames,1),2);
        [vec1, vec2] = anCueNames{idx,3};
        simil(ia) = cosSimil(vec1,vec2);
    end
killPermIU = sum(simil>conceptPermTH)/nperm; % what are the chances two images randomly chosen from our distribution of animals are at least as similar as the least similar concept neuron that we kicked out

permNumFP2_glob = zeros([1 numperm]);
x = zeros(1,10000);
x(1:killPermIU*10000) = 1;
for ih = 1:numperm % do this 10.000 times
    tic
    for numSU = 1:size(sensTrls_perm,2)
        perm = randi(numperm); % this actually might not be necessary
        
        ff = sum(ffIdx{numSU}(sensTrls_perm{numSU}(:, perm))); % take out a random permutation from that SU
        pp = sum(ppIdx{numSU}(sensTrls_perm{numSU}(:, perm)));
        
        %          && dpCor (rownum,perm) > threshDPcor(rownum) % threshDPcor would
        %          be a lower threshold because its not WC
        if ff > 0 && pp > 0
            %% chance that it is a concept neuron
            x = x(randperm(length(x)));
            if x(1) == 1
                continue
            end
            %%
            permNumFP2_glob(ih) = permNumFP2_glob(ih) + 1;
        end
    end
    toc
end
th_fp2_glob = prctile(permNumFP2_glob,95) % gives significance / th at 15
sumnum = sum(empNumFP2>=permNumFP2_glob) % gives 9996
empNumFP2 % at 21
cd(savingpath)
save('globalTest', 'th_fp2_glob', 'sumnum', 'empNumFP2', 'encVec_allSU', 'retVec_allSU', 'permNumFP2_glob')

end % end of function