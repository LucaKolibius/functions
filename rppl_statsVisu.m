%% ANALYSIS PART
clear
% load('X:\Luca\data\allSbj\allSpksHZ.mat', 'allSpks'); % I no longer do spike ripple coincidence
load('X:\Luca\data\allSbj\rpplBund_noSpkIntOrth.mat', 'rpplBund');
% numMW = 1; % for macro recording
numMW = 8; % for micro recording
for tw = 0 % :1:3 % only whole trial
    for maxMW = 1 % 0:1:1 % only wire with maximum value
        
        numRpplIdx = [];
        numRpplNdx = [];
        lenRpplNdx = [];
        lenRpplIdx = [];
        rpplDenIdx = [];
        rpplDenNdx = [];
%         coincIdx   = [];
%         coincNdx   = [];
        
        for bund = 1 : length(rpplBund)
            
            bidsID   = rpplBund(bund).bidsID;
            sesh     = rpplBund(bund).sesh;
            bundname = rpplBund(bund).bundname;
            idxTrl   = rpplBund(bund).idxTrlBund; % INDEXED TRIALS

%             % IF WE HAVE NO LOCATION INFORMATION, CONTINUE
%             if isempty(rpplBund(bund).location)
%                 continue
%             end
            
%             % IF BUNDLE HAS NO SPIKES OR IS NOT IN THE HIPPOCAMPUS, CONTINUE WITH
%             % THE NEXT ONE
%             if isempty(rpplBund(bund).idxTrlBund) || ~strcmp(rpplBund(bund).location, 'hipp')
%                 continue
%             end
            
%             idxSpk = and(and(strcmp({allSpks.bidsID}, bidsID), strcmp({allSpks.sesh}, sesh)), strcmp({allSpks.bundlename}, bundname));
%             idxSpk = find(idxSpk == 1);
                        
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
%             rpplSpkCoinc = zeros(trlNum,8);
            
            for mw = 1:numMW
                rppls = rpplBund(bund).rppls{mw}; % START + END OF THE RIPPLES ON THE CURRENT MICROWIRE
                    
                % lets see if this produces an error
%                 if isempty(rppls) % NOT SURE IF THIS IS OK
%                     continue
%                 end
                
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
                    
%                     for spk = 1 : size(idxSpk,2)
%                         
%                         spkTms    = allSpks(idxSpk(spk)).spks; % all spikes on the current bundle
%                         spkTrls   = spkTms(spkTms >= encTrig(trl,1) & spkTms <= encTrig(trl,2)) - encTrig(trl,1)+1; % the spikes that occur in this trial
%                         temp(spk) = mean(trlRpplDen(spkTrls)) / mean(trlRpplDen);
%                    
%                     end % END OF SPIKE LOOP 
                        % outside of ripple loop!
                    
                    lenRppl(trl,mw)      = mean(lenRpplAll);
                    rpplDen(trl,mw)      = mean(trlRpplDen);
%                     rpplSpkCoinc(trl,mw) = mean(temp);
                    
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
%                 rpplSpkCoinc = max(rpplSpkCoinc, [], 2);
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
            
%             %  INDEXED TRIALS - SPIKE RIPPLE COINCIDENCE
%             tmp = rpplSpkCoinc( idxTrl,:);
%             coincIdx = [coincIdx; tmp(:)];
            
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
            
%             %  NON-INDEXED TRIALS - SPIKE RIPPLE COINCIDENCE
%             tmp = rpplSpkCoinc(~idxTrl,:);
%             coincNdx = [coincNdx; tmp(:)];
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
        
%         % VISUALIZE RIPPLE COINCIDENCES
%         close all
%         figure('units','normalized','outerposition',[0 0 1 1]);
%                 
%         [p, ~] = ranksum(coincIdx, coincNdx, 'tail', 'right');
%         [p_perm] = perm_ranksum(coincIdx, coincNdx);
%         title({sprintf('Spike-Ripple Coincidences (p_{ranksum} = %.2f | p_(perm) = %.2f)', p, p_perm), ...
%             sprintf('Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(coincIdx), mean(coincIdx), median(coincIdx)), ...
%             sprintf('Non-Indexed Trials (#%d; M = %.2f; MED = %.2f)', length(coincNdx), mean(coincNdx), median(coincNdx))});
%         hold on
%         hand1 = histogram(coincIdx, 0:1:10, 'Normalization','probability');
%         hand2 = histogram(coincNdx, 0:1:10, 'Normalization','probability');
%         hand2.FaceColor = [1,0,0];
%         plot([mean(coincIdx) mean(coincIdx)], get(gca, 'YLim'), 'color', 'b', 'linew', 4, 'linestyle', '-');
%         plot([mean(coincNdx) mean(coincNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '-');
%         plot([median(coincIdx) median(coincIdx)], get(gca, 'YLim'), 'color', 'b',     'linew', 4, 'linestyle', '--');
%         plot([median(coincNdx) median(coincNdx)], get(gca, 'YLim'), 'color', [1 0 0], 'linew', 3, 'linestyle', '--');
%         legend('Indexed', 'Non-Indexed')
%         %         sgtitle(tits)
%         xlabel(tits);
        
        saveas(gcf,['X:\Luca\visu\rpplCoin_', saveTits], 'png' )
        close all

    end
end