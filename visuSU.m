% checkout the trialwindow in loadInDat
% have spksEnc and spksRet already inverted there
% SU 21 in the first index neuron session (indexNeuron = 1) is super tonic!
function visuSU
whereAmI;
global prePath;

savepath = [prePath, '\Luca\visu\SANDERISCONVINCED\'];
addpath(genpath([prePath, '\Luca\toolboxes\TREBER']));
load([prePath, 'Luca\data\allSbj\allSpksHZ_encTwo_retResp.mat'], 'allSpks')

for su =  [153 169 174 177 460] %[82 458]%1:size(allSpks,2)
    
    if su == 458
        continue
    end
    
    close all
    
    %% SETUP FIGURE
    mFigH = figure('units', 'normalized'); clf;
    %     set(gcf, 'units','normalized','outerposition',[0 0 1 1])
    
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
     
    if sum(allSpks(su).idxTrl) < 2
        continue % I'd like to have at least two reinstated trials
    end
    
    [hitsIdx, subjID, spksEnc, spksRet, sensTrlsFFPP, animalCues, WS, averageWS, ffTrls, ppTrls, ewp, associates] = loadInDat(allSpks, su);
    %     spksEnc = spksEnc'; spksRet = spksRet';
    
    % make spiketimes to number of spikes
    spksEncNum = cellfun(@length, ewp.enc) ./ ewp.twEnc;
    spksRetNum = cellfun(@length, ewp.ret) ./ ewp.twRet;
    spksEncNum = spksEncNum(hitsIdx); % get rid of misses
    spksRetNum = spksRetNum(hitsIdx); % get rid of misses
    
    %     [normENC, normRET]= normSpikeNumber_double(spksEncNum, spksRetNum);
    normENC = (spksEncNum - mean(spksEncNum)) / std(spksEncNum);
    normRET = (spksRetNum - mean(spksRetNum)) / std(spksRetNum);
    
    origENC = normENC; % original encoding  trial series (presentation order)
    origRET = normRET; % original retrieval trial series (matched to ENC presentation order)
    origDP  = origENC .* origRET; % original EWP
    
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
%     mFigH = figure('units', 'pixels');
    subplot(6,8,[1:6, 9:14, 17:22, 25:30 33:38])
    
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
    yhandRast = ylabel('Trial Number (matched)');
    ylim([0.5 counter1+0.5])
    yticks(0:10:counter1);
    xticks('');
    %     mAx = gca;
    %     mAx.YAxis.FontWeight = 'bold';
    %     mAx.XAxis.FontWeight = 'bold';
    %     mAx.FontSize = 16;
    set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',16,'FontWeight','Bold', 'LineWidth', 2);
    
    xlim([-2 15])
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
    line([0 0],[0 counter1+1], 'LineStyle','-', 'Color',[0.6 0.6 0.6 0.5], 'LineWidth',3)
    line([2 2],[0 counter1+1], 'LineStyle','-', 'Color',[0.6 0.6 0.6 0.5], 'LineWidth',3)
    
    line([0 0],[0 counter1+1], 'LineStyle','-', 'Color',[0.6 0.6 0.6 0.5], 'LineWidth',3)
    line([2 2],[0 counter1+1], 'LineStyle','-', 'Color',[0.6 0.6 0.6 0.5], 'LineWidth',3)
    
    %     for i = 1:size(sensTrlsFFPP,2)
    %         trialNum = sensTrlsFFPP(i);
    %         line([-1 5],[trialNum trialNum], 'LineStyle','-', 'Color',[0.6 0.6 0.6 0.5], 'LineWidth',2)
    %     end
    
    %% Trial Series
    subplot(6,8,[7:8, 15:16, 23:24, 31:32 39:40])
    hold on
    %                     plotting
    %     plot(origENC, 'linew',3, 'Color', mBlue, 'marker', 'o');
    %     plot(origRET, 'linew',3, 'Color', mOrange, 'marker', 'o');
    % plot(origDP, 'linew',3, 'Color', [0.47,0.67,0.19], 'marker', 'o');
    
    barHand = bar(origDP);
    barHand.FaceColor = (mBlue+mOrange)/2;
    
    %                     aesthetics
    xlim([0.5 length(hitsIdx)+0.5])
    curY = ylim; 
    ylim([curY(1) curY(2) + 1]);
    %         hand = legend({'ENC', 'RET', ''});
    ylabel('Reinstatement value');
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
    set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',16,'FontWeight','Bold', 'LineWidth', 2);
    box on
    
    %     % this marks only the indexed trials
    %     recPos = get(gca,'Ylim');
    %     startP = -0.5;
    %     for iy = 1 : size(normENC,2)
    %         startP = startP+1;
    %         if any(1 == allSpks(su).idxTrl)
    %             r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[0.6 0.6 0.6 0.5]); % dark
    %             r.EdgeColor = [1 1 1]; % make edgecolor white
    %             %                 elseif ppIdx{suNum}(iy)==1 && sensTrls(iy) ==1
    %             %                     r = rectangle('Position',[startP recPos(1) 1 abs(recPos(1))+abs(recPos(2))],'FaceColor',[1 1 1 0.35]); % white
    %         end
    %     end
    
    
    %% fplot
    subplot(6,8,[41:46])
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
    plot(ploDT, encConv, 'linew', 3, 'color', mBlue);
    plot(ploDT, retConv, 'linew', 3, 'color', mOrange);
    %     ylabel({'Smoothed', 'Occurences', ''})
    
    box on
    xlabel('Time after cue onset (seconds)')
    xticks([-2 0 2 5:5:20]);
    yticks('')
    xlim([-2 15])
    
    yhandFire = ylabel({'Spike density'});
    newPos    = yhandFire.Position;
    newPos(1) = yhandRast.Position(1); % align on x axis with rasterplot ylabel
    yhandFire.Position = newPos;
      
    xticklabels({'-2', 'Cue Onset', '2', '5', '10', '15', '20'})
    
    %     mAx = gca;
    %     mAx.YAxis.FontWeight = 'bold';
    %     mAx.XAxis.FontWeight = 'bold';
    %     mAx.FontSize = 14;
    set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',16,'FontWeight','Bold', 'LineWidth', 2);
    title('Firing rate during reinstated trials', 'FontName', 'Arial', 'FontWeight', 'Bold', 'FontSize', 14)

    line([0 0],[get(gca, 'YLim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3) % cue onset
    line([2 2],[get(gca, 'YLim')], 'LineStyle','-', 'Color',[0.6 0.6 0.6], 'LineWidth',3) % stimulus onset
    
%     hold on
%     L(1) = plot(nan, nan, 'color', mBlue);
%     L(2) = plot(nan, nan, 'color', mOrange);
%     [legPos, hobj, ~, ~] = legend(L, {'Encoding', 'Retrieval'}, 'FontSize',16, 'FontWeight', 'bold');
%     hl = findobj(hobj,'type','line');
%     set(hl,'LineWidth',15);
%     set(legPos, 'Position', [0.033 0.875 0.078 0.061])
%     legend('boxoff');
    
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
    %     figure
    %     subplot(6,8,[39:40])
    subplot(6,8,[47 48])
    
    density_plot(WS, [1:size(WS,1)], 1);

    %%
    subplot(6,8,[41:46])
    yhandFire = ylabel({'Spike density'});
    newPos = yhandFire.Position;
    newPos(1) = yhandRast.Position(1); % align on x axis with rasterplot ylabel
    yhandFire.Position = newPos;

    %% animal cue
    %  get cosine similarity between two animals
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
    
    
    load([prePath, 'Luca\old stuff\anCueNames.mat'], 'anCueNames');
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
    
    chooseAni = animalCues';
    chooseAni = [chooseAni, associates];
    chooseAni = chooseAni';
    
    %% SAVING
%     saveas(gcf, [savepath, 'visuIU_',subjID, '_SU',num2str(su)], 'svg'); % emf, svg
        
%% ANIMAL CUE AND ASSOCIATE IMAGES
    figure(2); clf;
    cntr = 1;
    % UNCOMMENT IF YOU WANT TO HAVE THE ANIMAL CUE
    for anCue = 1:numel(chooseAni)
%         if anCue > 2
%             continue
%         end
        
%         anCuePos = [47, 48];
                subplot(size(chooseAni,2),3,cntr)
        
%         subplot(6,8,anCuePos(anCue))
        
        ylim([0 1])
        xlim([0 1])
        imName = char(chooseAni(anCue));
        
        imHandle = imread(imName);
        image(imHandle, 'XData', [-0.1+0.1*anCue 0.1*anCue], 'YData', [0 1]);
        box off
        
        axis off
        axis square
        
        cntr = cntr + 1;
    end
    
    %     title(sprintf('Similarity: %.3f (90%%-th: %.3f)', minVal, word2vec_th))
    %     mAx = gca;
    %     mAx.Title.Position = [0.0500 1.5 0];
    
    %% saving image
    cd(savepath)
    
%     % maximize figure on screen. adapted for 1 to N monitors.
%     %    % Get pixel position of monitors.
%     MP = get(0, 'MonitorPositions');
%     N = size(MP, 1);
%     % Might want to set an initial position this to some reasonable location
%     % in the event of the window being "Restored Down".
%     newPosition = MP(1,:);
%     
%     if size(MP, 1) == 1
%         % Single monitor
%         set(mFigH,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
%     else
%         % Multiple monitors - shift to the Nth monitor.
%         newPosition(1) = newPosition(1) + MP(N,1);
%     end
%     mFigH.set('Position', newPosition, 'units', 'normalized');
%     mFigH.WindowState = 'maximized'; % Maximize with respect to current monitor.
    
    saveas(gcf, [savepath, 'visuIUcue_',subjID, '_SU',num2str(su)], 'png'); % emf, svg
    %     catch err
    %         fprintf(1,'The identifier was:\n%s',err.identifier);
    %         fprintf(2,'There was an error! The message was:\n%s',err.message);
    %         disp('Catching.');
    %         continue
    %     end
end
    cd(savepath)

end % end of function

% tt = [allSpks.iu];
% tt(tt==1) = 0;
% find(tt==2)