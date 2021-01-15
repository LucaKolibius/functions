% checkout the trialwindow in loadInDat
% have spksEnc and spksRet already inverted there
% SU 21 in the first index neuron session (indexNeuron = 1) is super tonic!
function visuSU
whereAmI(0);
global prePath;

savepath = [prePath, '\Luca\SFN2020'];
addpath(genpath([prePath, '\Luca\toolboxes\TREBER']));
load([prePath, 'Luca\data\allSbj\allSpksHZ.mat'], 'allSpks')

for su = 405%[177, 405]%1:size(allSpks,2)
    close all



    if sum(allSpks(su).idxTrl) == 0
        continue % NO INDEX TRIALS
    end
    
    [hitsIdx, subjID, spksEnc, spksRet, sensTrlsFFPP, animalCues, WS, averageWS, ffTrls, ppTrls] = loadInDat(allSpks, su);
    spksEnc = spksEnc'; spksRet = spksRet';
    
    % make spiketimes to number of spikes
    spksEncNum = num2cell(cellfun(@length, spksEnc));
    spksRetNum = num2cell(cellfun(@length, spksRet));
    spksEncNum = spksEncNum(:,hitsIdx); % get rid of misses
    spksRetNum = spksRetNum(:,hitsIdx); % get rid of misses
    
    [~ , normENC, normRET] = normSpikeNumber(spksEncNum, spksRetNum);
    
    origENC = normENC; % original encoding  trial series (presentation order)
    origRET = normRET; % original retrieval trial series (presentation order)
    origDP  = origENC .* origRET; % original dot-product
    
    %% SCALE RASTERPLOT BY EWP
    maxDP = max(origDP);
    scalDP = origDP/maxDP;
    scalDP(scalDP<0) = 0;
    scalDP = scalDP + 0.2;
    scalDP(scalDP>1) = 1;
    
    %% DON'T SCALE
%     scalDP = ones(size(origENC));
    
    %% Visualization
    timeWindow = [-3 30];
    dt = linspace(timeWindow(1), timeWindow(2), (abs(timeWindow(1)) + abs(timeWindow(2))) *1000+1 );
    ploDT = dt+0.0005; ploDT(end) = [];
    % rasterplot
    mFigH = figure('units', 'pixels');
    subplot(6,8,[1:6, 9:14, 17:22, 25:30])
    
%     mainTitle = sprintf(['Participant: ', subjID, 'Session: ', allSpks(su).sesh, '  |  SU#%i'], allSpks(su).su);
%     sgtitle(mainTitle); % we need at least R2018b
    
    hold on
    n_encHit = [];
    counter1 = 0;
    
    % % Encoding
%     mBlue = [0 0.4470 0.7410];
mBlue = [[0.368627450980392,0.235294117647059,0.600000000000000]];

    for trl = 1:size(hitsIdx,1) % trials
        counter1 = counter1+1;
        
        x = spksEnc{1, hitsIdx(trl)};
        xd = [x';x'];
        y = counter1*ones(1,length(x));
        y = [y-.5;y+.5];
        if ~isempty(x) % so I don't clear the handle, which will create problems with the legend later
            %             lineH{1} = line(xd,y,'Color',[mBlue, scalDP(trl),'LineWidth',2);
            plot(xd,y,'Color',[mBlue, scalDP(trl)],'LineWidth',2);
        end
        [n_encHit(counter1,:),~] = histcounts(x,dt);
    end
    
    %     title('Rasterplot encoding(blue) + retrieval(red)');
    ylabel('Trial Number (matched)');
    ylim([0.5 counter1+0.5])
    yticks(0:10:counter1);
    xticks('');
    mAx = gca;
    mAx.YAxis.FontWeight = 'bold';
    mAx.XAxis.FontWeight = 'bold';
    mAx.FontSize = 16;
    xlim([-2 20])
    box on
    
    % % RETRIEVAL
%     mOrange = [0.8500 0.3250 0.0980];
mOrange = [[0.901960784313726,0.380392156862745,0.00392156862745098]];

    hold on
    counter1 = 0;
    n_retHit = [];
    for trl = 1 : size(hitsIdx,1) % trials
        counter1 = counter1+1;
        x = spksRet{1, hitsIdx(trl)};
        xd = [x';x'];
        y = counter1*ones(1,length(x));
        y = [y-.5;y+.5];
        if ~isempty(x)
            %             lineH{3} = line(xd,y,'Color', mOrange,'LineWidth',2);
            plot(xd,y,'Color',[mOrange, scalDP(trl)],'LineWidth',2);
        end
        [n_retHit(counter1,:),~] = histcounts(x,dt);
    end
    line([0 0],[0 counter1+1], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3)
    line([2 2],[0 counter1+1], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3)
    
    %     for i = 1:size(sensTrlsFFPP,2)
    %         trialNum = sensTrlsFFPP(i);
    %         line([-1 5],[trialNum trialNum], 'LineStyle','-', 'Color',[0.6 0.6 0.6 0.5], 'LineWidth',2)
    %     end
    
    %% Trial Series
    subplot(6,8,[7:8, 15:16, 23:24, 31:32])
    hold on
    %                     plotting
%     plot(origENC, 'linew',3, 'Color', mBlue, 'marker', 'o');
%     plot(origRET, 'linew',3, 'Color', mOrange, 'marker', 'o');
% plot(origDP, 'linew',3, 'Color', [0.47,0.67,0.19], 'marker', 'o');

barHand = bar(origDP);
barHand.FaceColor = (mBlue+mOrange)/2;

    %                     aesthetics
    xlim([0.5 length(hitsIdx)+0.5])
    %         hand = legend({'ENC', 'RET', ''});
    ylabel('Standardized firing rate');
    xticks('')
    %         xlabel('Trial Number')
    %             title([subjID, sprintf('-SU%.0f', su)]);
    %         set(mFigH,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
    axH = gca;
    axH.YAxis.FontWeight = 'bold';
    axH.XAxis.FontWeight = 'bold';
    axH.FontSize = 16;
    axH.XDir = 'reverse';
    axH.XAxisLocation = 'top';
    axH.YAxisLocation = 'right';
    view([90 90])
    
    % this marks only the indexed trials
    recPos = get(gca,'Ylim');
    startP = -0.5;
    for iy = 1 : size(normENC,2)
        startP = startP+1;
        if any(iy == allSpks(su).idxTrl)
            r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[0.6 0.6 0.6 0.5]); % dark
            r.EdgeColor = [1 1 1]; % make edgecolor white
            %                 elseif ppIdx{suNum}(iy)==1 && sensTrls(iy) ==1
            %                     r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[1 1 1 0.35]); % white
        end
    end
    
    
    %% fplot
    subplot(6,8,[33:38, 41:46])
    hold on
    
    retIdxSpikes = n_retHit(allSpks(su).idxTrl,:);
    retIdxSpikes = sum(retIdxSpikes,1);
    encIdxSpikes = n_encHit(allSpks(su).idxTrl,:);
    encIdxSpikes = sum(encIdxSpikes,1);
    
    % gaussian kernel
    mlength = [-0.075:0.0002:0.075];
    mSigma  = 0.02; % 20ms
    % mSigma  = 0.01;
    mKernel = normpdf(mlength,0,mSigma);
    mKernel = mKernel/max(mKernel); % normalize peak to 1
    
    % encoding
    encConv   = conv(mKernel, encIdxSpikes);
    % get rid of edges
    encConv(1:375)       = [];
    encConv(end-374:end) = [];
    
    % retrieval
    retConv = conv(mKernel,retIdxSpikes);
    % get rid of edges
    retConv(1:375)       = [];
    retConv(end-374:end) = [];
    
    % plot convolved series
    hold on
    title('Firing rate during indexed trials')
    plot(ploDT, encConv, 'linew', 3, 'color', mBlue);
    plot(ploDT, retConv, 'linew', 3, 'color', mOrange);
    ylabel({'Smoothed', 'Occurences', ''})
    box on
    xlabel('Time after cue onset (seconds)')
    xticks([-2 0 2 5:5:20]);
    yticks('')
    xlim([-2 20])

    xticklabels({'-2', 'Cue Onset', '2', '5', '10', '15', '20'})
    
    mAx = gca;
    mAx.YAxis.FontWeight = 'bold';
    mAx.XAxis.FontWeight = 'bold';
    mAx.FontSize = 14;
    line([0 0],[get(gca, 'YLim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3) % cue onset
    line([2 2],[get(gca, 'YLim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3) % stimulus onset
    
    hold on
    L(1) = plot(nan, nan, 'color', mBlue);
    L(2) = plot(nan, nan, 'color', mOrange);
    [legPos, hobj, ~, ~] = legend(L, {'Encoding', 'Retrieval'}, 'FontSize',16, 'FontWeight', 'bold');
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',15);
    set(legPos, 'Position', [0.033 0.875 0.078 0.061])
    legend('boxoff');

    %     %% frequency plot
    %     sp2H = subplot(4,1,2);
    %     title('Frequency Plot')
    %     hold on
    %     frH = sum(n_encHit,1)./size(n_encHit,1)./0.125; % I have a bin each 250ms/transforms into herz
    %     frH([1 end]) = []; % cut off the wings
    %     dtEnc([1 end]) = [];
    %     plot(dtEnc,frH,'ks-','LineWidth',3, 'Color', mBlue);
    %     xlim([-1 5])
    %     [legPos, ~, ~, ~] = legend({'Enc - Hits', 'Enc - Misses'}, 'FontSize',8, 'FontWeight','bold');
    %     set(legPos, 'Position', [0.85 0.6 0.06 0.03])
    %     legend('boxoff');
    %     xlabel('Time in Seconds')
    %     ylabel('Frequency [Hz]');
    %     mAxA = gca;
    %     mAxA.YAxis.FontWeight = 'bold';
    %     mAxA.XAxis.FontWeight = 'bold';
    %     mAxA.FontSize = 12;
    %     mAxA.XAxis.FontWeight = 'bold';
    %     curYlimA = get(mAxA, 'YLim');
    %
    %     %
    %     frH = sum(n_retHit,1)./size(n_retHit,1)./0.125; % I have a bin each 250ms/transforms into herz
    %     frH([1 end]) = []; % cut off the wings
    %     dtRet([1 end]) = []; % cut off the wings
    %     plot(dtRet,frH,'ks-','LineWidth',3, 'Color', mOrange);
    %     [legPos, ~, ~, ~] = legend({'Ret - Hits', 'Ret - Misses'}, 'FontSize',8, 'FontWeight','bold');
    %     set(legPos, 'Position', [0.8462 0.6326 0.0594 0.0337])
    %     legend('boxoff');
    %
    %
    %     xlabel('Time in Seconds')
    %     ylabel('Frequency [Hz]');
    %     mAxB = gca;
    %     mAxB.YAxis.FontWeight = 'bold';
    %     mAxB.XAxis.FontWeight = 'bold';
    %     mAxB.FontSize = 12;
    %     mAxB.XAxis.FontWeight = 'bold';
    %     curYlimB = get(mAxB, 'YLim'); % sets minimum to 0 % commenting this out leads to the retrieval frequency plot using the same y-axis as the encoding freqp
    %
    %     % set the YLim value for the both frequency plots to the maximum of either. This will lead to equal YLims as well as no cut-offs
    %     mAxA.YLim = [0 max(curYlimA(2), curYlimB(2)) ];
    %     mAxB.YLim = [0 max(curYlimA(2), curYlimB(2)) ];
    %
    %     % draw grey lines
    %     line([0 0],[get(sp2H,'YLim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3)
    %

    %% waveshape
    figure
%     subplot(6,8,[39:40])
    density_plot(WS, [1:size(WS,1)], 1);
    title(sprintf('Waveshapes (%d Spikes)', size(WS,1)))
    xlabel('Seconds')
    caYticks = yticks;
    try
    yticks([caYticks(1) 0 caYticks(end)]) % this produces an error if the first tick in the current yticks is a 0 already making this [0 0 X];
    catch
    end
    
    axH = gca;
    axH.YAxis.FontWeight = 'bold';
    axH.XAxis.FontWeight = 'bold';
    axH.FontSize = 12;
    axH.YLabel.VerticalAlignment = 'middle';
    axH.XLabel.VerticalAlignment = 'cap';
    
    %% animal cue
%     load([prePath, '\Luca\data\allSbj\anCueNames.mat'], 'anCueNames');
%     nperm = 1;
%     
%     simil = zeros(1,nperm);
%     for ia=1:nperm
%         idx = randperm(size(anCueNames,1),2);
%         [vec1, vec2] = anCueNames{idx,3};
%         simil(ia) = cosSimil(vec1,vec2);
%     end
%     word2vec_th = prctile(simil,word2vec_plvl);
%     
    %
%     load('X:\Luca\anCueNames.mat', 'anCueNames');
    numSignTrls = size(animalCues,2);
    if strcmp(animalCues{1}(end-2:end),'bmp')
        %         cd('Z:\Luca\Github\FacePlace-task\EMpairs_v5_2017-09-11\image_data\EMtune\fVSp_resized\empairs_cropped');
        cd([prePath, '\Luca\exp_code\FacePlace-task\EMpairs_v5_2017-09-11\image_data\EMtune\fVSp_resized\empairs_cropped']);
    elseif strcmp(animalCues{1}(end-2:end),'jpg')
        %         cd('Z:\Luca\Github\FacePlace-task\EMpairs_v5_2017-09-11\image_data\EMtune\formatted_180-180')
%         cd('X:\Luca\expStimuli\formatted_180-180')
        cd([prePath, '\Luca\exp_code\FacePlace-task\EMpairs_v5_2017-09-11\image_data\EMtune\formatted_180-180']);
    end
    
%     % find the most different animal cues if there are more than two trials
%     empSimi = [];
%     comp = [];
%     for it = 1:numSignTrls
%         vec1 = anCueNames{ismember(anCueNames(:,1), animalCues(it)),3};
%         for iz = 1:numSignTrls
%             if iz == it || iz<it
%                 continue
%             end
%             vec2 = anCueNames{ismember(anCueNames(:,1), animalCues(iz)),3};
%             empSimi = [empSimi cosSimil(vec1, vec2)];
%             
%             comp = [comp, [{animalCues(it)}; {animalCues(iz)}]];
%         end
%     end
%     [minVal,minIdx] = min(empSimi); % the two animals with the greatest distance
%     chooseAni = [comp{:,minIdx}];
    
chooseAni = animalCues;

% UNCOMMENT IF YOU WANT TO HAVE THE ANIMAL CUE
    for anCue = 1:size(chooseAni,2)
        if anCue > 2
            continue
        end
        
        anCuePos = [47, 48];
        subplot(6, 8, anCuePos(anCue))
        ylim([0 1])
        xlim([0 1])
        imName = char(chooseAni(anCue));
        
        imHandle = imread(imName);
        image(imHandle, 'XData', [-0.1+0.1*anCue 0.1*anCue], 'YData', [0 1]);
        box off
        
        axis off
        axis square
    end

%     title(sprintf('Similarity: %.3f (90%%-th: %.3f)', minVal, word2vec_th))
%     mAx = gca;
%     mAx.Title.Position = [0.0500 1.5 0];
    
    %% saving image
    cd(savepath)
    
    % maximize figure on screen. adapted for 1 to N monitors.
    %    % Get pixel position of monitors.
    MP = get(0, 'MonitorPositions');
    N = size(MP, 1);
    % Might want to set an initial position this to some reasonable location
    % in the event of the window being "Restored Down".
    newPosition = MP(1,:);
    
    if size(MP, 1) == 1
        % Single monitor
        set(mFigH,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
    else
        % Multiple monitors - shift to the Nth monitor.
        newPosition(1) = newPosition(1) + MP(N,1);
    end
    mFigH.set('Position', newPosition, 'units', 'normalized');
    mFigH.WindowState = 'maximized'; % Maximize with respect to current monitor.
    
    saveas(gcf, ['BvisuSU_',subjID, '_SU',num2str(su)], 'svg'); % emf
%     catch err
%         fprintf(1,'The identifier was:\n%s',err.identifier);
%         fprintf(2,'There was an error! The message was:\n%s',err.message);
%         disp('Catching.');
%         continue
%     end
end
end % end of function

% tt = [allSpks.iu];
% tt(tt==1) = 0;
% find(tt==2)