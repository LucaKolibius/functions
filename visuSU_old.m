% adapt rasterplot linewidth based on number of spikes
% density waveshape

fetch_allSubj
for sbj = 1 : size(allSubj, 1)
    subjID = allSubj{sbj};
    
    
    % cd
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
    abc = dir; cd(abc(3).name);
    
    
    % logfile
    p2d = cd;
    p2d(end+1)='\';
    [~, hitsIdx, missIdx, ~, ~, retTrigger, encTrigger] = loadLogs(p2d);
    
    % load single units
    cd advancedAnalysis\postXCrej
    loadVar = dir('tableTimestamps_postXC_*'); load(loadVar.name); % load tableTimestamps_postXC
    cd ../manualRej/
    loadVar = dir('clNames*'); load(loadVar.name); % load clNames
    cd ../..
    
    clCount = 0;
    allSpikes = [];
    for su=1:size(tableTimestamps_postXC,2)
        mVarname = tableTimestamps_postXC.Properties.VariableNames{1,su};
        if ismember(mVarname, clNames)
            clCount = clCount +1;
            clusterSpikes = table2array(tableTimestamps_postXC{1,su})'; % load in all spiketimes from cluster mVarname
            clusterSpikes = clusterSpikes/32000; % samples to sec
            
            % segment trials of this SU into trials
            spksEnc(:, clCount) = insertSpiketimes(encTrigger, clusterSpikes, 1, [-1 5])'; % cue locked; 6s  %I could also generate a big file with all subj data in it (using a struct) or put enc and ret in one variable, but what for?
            spksRet(:, clCount) = insertSpiketimes(retTrigger, clusterSpikes, 1, [-1 3])'; % cue locked; 4s
            
            % get waveshape
            demGate = cd; % where the ncs files are
            
            if ~isempty(regexp(mVarname, 'Pos', 'ONCE'))
                posNeg = 'posDetect';
            elseif ~isempty(regexp(mVarname, 'Neg', 'ONCE'))
                posNeg = 'negDetect';
            end
            
            cd(posNeg)
            
            temp1 = mVarname;
            temp1(end-5:end) = [];
            temp = dir(['times_CSC_', temp1, '.mat']);
            
            load(temp.name, 'spikes', 'cluster_class');
            
            clNum = mVarname(end);
            extrCl = cluster_class(:,1)==str2num(clNum);
            
            counterX = 0;
            for ib = 1 : size(extrCl,1)
                if extrCl(ib,1) == 1
                    counterX = counterX+1;
                    allSpikes(clCount).waveshape(counterX, :) = spikes(ib,:);
                end
            end
            averageWS(clCount,:) = mean(allSpikes(clCount).waveshape,1);
            cd(demGate)
        end
    end
    
    %% Visualization
    timeWindow = [-1 5];
    dtEnc = linspace(timeWindow(1),timeWindow(2),49); % steps of 125
    timeWindow = [-1 3];
    dtRet = linspace(timeWindow(1),timeWindow(2), (abs(timeWindow(1)) + abs(timeWindow(2))) *8+1 ); % steps of 125ms
    
    for su = 1 : size(spksEnc,2) % single units
        % rasterplot
        clf % clear current figure window
        mFigH = figure(93);
        subplot(3,2,1)
        
        session = subjID;
        session(regexp(subjID, '_'))='-'; % switch _ to -
        mainTitle = sprintf(['Session: ', session, '  |  SU#%i'], su);
        %         sgtitle(mainTitle); % we need at least R2018b
        
        hold on
        n_encHit = [];
        counter1 = 0;
        
        % % HITS
        for trl = 1:size(hitsIdx,1) % trials
            counter1 = counter1+1;
            x = spksEnc{hitsIdx(trl),su};
            xd = [x;x];
            y = counter1*ones(1,length(x));
            y = [y-.5;y+.5];
            if ~isempty(x) % so I don't clear the handle, which will create problems with the legend later
                lineH{1} = line(xd,y,'Color','b','LineWidth',2);
            end
            [n_encHit(hitsIdx(trl),:),~] = hist(x,dtEnc);
        end
        line([-1 5], [counter1 counter1], 'LineStyle', '--', 'Color','k', 'LineWidth',1)
        
        % % MISS
        n_encMiss = [];
        for trl = 1 : size(missIdx,1) % trials
            counter1 = counter1+1;
            x = spksEnc{missIdx(trl),su};
            xd = [x;x];
            y = counter1*ones(1,length(x));
            y = [y-.5;y+.5];
            if ~isempty(x)
                lineH{2} = line(xd,y,'Color','r','LineWidth',2);
            end
            [n_encMiss(missIdx(trl),:),~] = hist(x,dtEnc);
        end
        line([0 0],[0 counter1], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3)
        title('Encoding');
        ylabel('Trial Number (re-sorted)');
        ylim([0 counter1+1])
        yticks(0:20:counter1+1); % this makes problems
        mAx = gca;
        mAx.YAxis.FontWeight = 'bold';
        mAx.XAxis.FontWeight = 'bold';
        mAx.FontSize = 12;
        xlim([-1 5])
%         box on
        
        % % RETRIEVAL - HIT
        subplot(3,2,2)
        hold on
        counter1 = 0;
        n_retHit = [];
        for trl = 1 : size(hitsIdx,1) % trials
            counter1 = counter1+1;
            x = spksRet{hitsIdx(trl),su};
            xd = [x;x];
            y = counter1*ones(1,length(x));
            y = [y-.5;y+.5];
            if ~isempty(x)
                lineH{3} = line(xd,y,'Color','c','LineWidth',2);
            end
            [n_retHit(hitsIdx(trl),:),~] = hist(x,dtRet);
        end
        line([-1 3], [counter1 counter1], 'LineStyle', '--', 'Color','k', 'LineWidt',1)
        
        % % RETRIEVAL - MISS
        n_retMiss = [];
        for trl = 1 : size(missIdx,1) % trials
            counter1 = counter1+1;
            x = spksRet{missIdx(trl),su};
            xd = [x;x];
            y = counter1*ones(1,length(x));
            y = [y-.5;y+.5];
            if ~isempty(x)
                lineH{4} = line(xd,y,'Color','m','LineWidth',2);
            end
            [n_retMiss(missIdx(trl),:),~] = hist(x,dtRet);
        end
        line([0 0],[0 counter1], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3)
        
        ylabel('Trial Number (re-sorted)');
        ylim([0 counter1+1])
        yticks(0:20:counter1+1); % this makes problems
        
        mAx = gca;
        mAx.YAxis.FontWeight = 'bold';
        mAx.XAxis.FontWeight = 'bold';
        mAx.FontSize = 12;
        xlim([-1 3])
        title('Retrieval');
%         box on
        
        %         hold on
        %         L(1) = plot(nan, nan, 'b');
        %         L(2) = plot(nan, nan, 'r');
        %         L(3) = plot(nan, nan, 'c');
        %         L(4) = plot(nan, nan, 'm');
        %         [legPos, hobj, ~, ~] = legend(L, {'Enc - Hits', 'Enc - Misses', 'Ret - Hits', 'Ret - Misses'}, 'FontSize',8, 'FontWeight','bold');
        %         hl = findobj(hobj,'type','line');
        %         set(hl,'LineWidth',2);
        %         set(legPos, 'Position', [0.91, 0.86, 0.06, 0.064])
        %         legend('boxoff');
        
        %% frequency plot
        sp2H = subplot(3,2,3);
        hold on
        frH = sum(n_encHit,1)./size(n_encHit,1)./0.125; % I have a bin each 250ms/transforms into herz
        frH([1 end]) = []; % cut off the wings
        frM = sum(n_encMiss,1)./size(n_encMiss,1)./0.125;
        frM([1 end]) = []; % cut off the wings
        dtEnc([1 end]) = [];
        plot(dtEnc,frH,'ks-','LineWidth',3, 'Color', 'b');
        plot(dtEnc,frM,'ks-','LineWidth',3, 'Color', 'r');
        xlim([-1 5])
        [legPos, ~, ~, ~] = legend({'Enc - Hits', 'Enc - Misses'}, 'FontSize',8, 'FontWeight','bold');
        set(legPos, 'Position', [0.4050 0.6326 0.0594 0.0337])
        legend('boxoff');
        xlabel('Time in Seconds')
        ylabel('Frequency [Hz]');
        mAxA = gca;
        mAxA.YAxis.FontWeight = 'bold';
        mAxA.XAxis.FontWeight = 'bold';
        mAxA.FontSize = 12;
        mAxA.XAxis.FontWeight = 'bold';
        curYlimA = get(mAxA, 'YLim');
        
        sp2bH = subplot(3,2,4);
        hold on
        frH = sum(n_retHit,1)./size(n_retHit,1)./0.125; % I have a bin each 250ms/transforms into herz
        frH([1 end]) = []; % cut off the wings
        frM = sum(n_retMiss,1)./size(n_retMiss,1)./0.125;
        frM([1 end]) = []; % cut off the wings
        dtRet([1 end]) = []; % cut off the wings
        plot(dtRet,frH,'ks-','LineWidth',3, 'Color', 'c');
        plot(dtRet,frM,'ks-','LineWidth',3, 'Color', 'm');
        [legPos, ~, ~, ~] = legend({'Ret - Hits', 'Ret - Misses'}, 'FontSize',8, 'FontWeight','bold');
        set(legPos, 'Position', [0.8462 0.6326 0.0594 0.0337])
        legend('boxoff');
        
        
        xlabel('Time in Seconds')
        ylabel('Frequency [Hz]');
        mAxB = gca;
        mAxB.YAxis.FontWeight = 'bold';
        mAxB.XAxis.FontWeight = 'bold';
        mAxB.FontSize = 12;
        mAxB.XAxis.FontWeight = 'bold';
        curYlimB = get(mAxB, 'YLim'); % sets minimum to 0 % commenting this out leads to the retrieval frequency plot using the same y-axis as the encoding freqp
        
        % set the YLim value for the both frequency plots to the maximum of either. This will lead to equal YLims as well as no cut-offs
        mAxA.YLim = [0 max(curYlimA(2), curYlimB(2)) ];
        mAxB.YLim = [0 max(curYlimA(2), curYlimB(2)) ];
        
        % draw grey lines
        sp2H = subplot(3,2,3);
        line([0 0],[get(sp2H,'YLim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3)
        
        sp2bH = subplot(3,2,4);
        line([0 0],[get(sp2bH,'YLim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3)
        
        %% bar plot for average frequency
        % encoding
        sp3H = subplot(325);
        hold on
        numEnc = cellfun(@length,spksEnc(:,su)); % number of spikes per trial
        numEnc_hits = numEnc(hitsIdx)/7;
        numEnc_miss = numEnc(missIdx)/7;
        
        numEnc_hitsMean  = mean(numEnc_hits);
        numEnc_hitsStd   = std(numEnc_hits);
        
        numEnc_missMean  = mean(numEnc_miss);
        numEnc_missStd   = std(numEnc_miss);
        
        hold on
        bar(1,numEnc_hitsMean, 'Blue')
        errorbar(1,numEnc_hitsMean, numEnc_hitsStd, numEnc_hitsStd, 'Color', 'k', 'LineWidth', 1.5);
        bar(2,numEnc_missMean, 'Red')
        errorbar(2,numEnc_missMean, numEnc_missStd, numEnc_missStd, 'Color', 'k', 'LineWidth', 1.5);
        
        %% repeat for retrieval
        numRet = cellfun(@length,spksRet(:,su)); % number of spikes per trial
        numRet_hits = numRet(hitsIdx)/5;
        numRet_miss = numRet(missIdx)/5;
        
        numRet_hitsMean  = mean(numRet_hits);
        numRet_hitsStd   = std(numRet_hits);
        
        numRet_missMean  = mean(numRet_miss);
        numRet_missStd   = std(numRet_miss);
        
        hold on
        bar(3,numRet_hitsMean, 'Cyan')
        errorbar(3,numRet_hitsMean, numRet_hitsStd, numRet_hitsStd, 'Color', 'k', 'LineWidth', 1.5);
        bar(4,numRet_missMean, 'Magenta')
        errorbar(4,numRet_missMean, numRet_missStd, numRet_missStd, 'Color', 'k', 'LineWidth', 1.5);
        
        %
        xlim([0.5 4.5])
        xticks(1:4)
        xticklabels({'Encoding - Hit', 'Encoding - Miss', 'Retrieval - Hit', 'Retrieval - Miss'})
        ylabel('Firing Rate [Hz]', 'FontSize', 8, 'FontWeight', 'bold');
        sp3H.YGrid = 'On';
        sp3H.YMinorGrid = 'On';
        box off
        mAx = gca;
        mAx.YAxis.FontWeight = 'bold';
        mAx.XAxis.FontWeight = 'bold';
        mAx.FontSize = 12;
        
        %% waveshape
        subplot(326)
        hold on
        
        for ib = 1 : size(allSpikes(su).waveshape,1)
            plot(allSpikes(su).waveshape(ib, :),'k')
        end
        plot(averageWS(su,:), 'r');
        
        
        xlim([1 64])
        ylabel('Amplitude in \muVolt')
        xticks([]);
        mAx = gca;
        mAx.YAxis.FontWeight = 'bold';
        mAx.XAxis.FontWeight = 'bold';
        mAx.FontSize = 12;
        
        %% saving image
        cd advancedAnalysis/
        mkdir('visuSU'); cd visuSU;
        set(mFigH,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
        saveas(gcf, ['visuSU_',subjID, '_SU',num2str(su)], 'jpg');
        
        
    end
end



