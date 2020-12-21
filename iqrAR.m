function [isAR, artSnip] = iqrAR(LFP, mkPlot)

numChan = size(LFP,1);       % number of channel
lenLFP  = size(LFP,2);       % number of time point in LFP
isAR    = zeros(numChan,lenLFP);  % preallocate output
gradLFP = diff(LFP,1,2);     % difference LFP

for chan = 1 : numChan
    
    chanLFP = LFP(chan,:);
    chanGra = gradLFP(chan,:);
    thH = median(chanLFP) + 3*iqr(chanLFP); % iqr = prctile(x,75) - prctile(x,25)
    thL = median(chanLFP) - 3*iqr(chanLFP); % iqr = prctile(x,75) - prctile(x,25)
    
    
    %% ARTEFACTS OVER THE THRESHOLD
    ovTh       = or(chanLFP>thH, [0 chanGra>thH]);
    ovTh_diff  = diff(ovTh);
    ovTh_start = find(ovTh_diff ==  1)'; % start of the artefact
    ovTh_end   = find(ovTh_diff == -1)'; % end of the artefact
    
    
    
%     ovTh       = chanGra>thH;
%     ovTh_diff  = diff(ovTh);
%     ovTh_start = [ovTh_start; find(ovTh_diff ==  1)'];
%     ovTh_end   = [ovTh_end  ; find(ovTh_diff == -1)'];
    
    % there will be a problem when an artefact precedes and exceeds
    % the snippet because the start and end lengths will be the same
    if length(ovTh_end) == length(ovTh_start)-1 % artefact extends the snippet
        ovTh_end = [ovTh_end; lenLFP];
    end
    
    if length(ovTh_end)-1 == length(ovTh_start) % artefact precedes the snippet
        ovTh_start = [1; ovTh_start];
    end
    
    ovThArt    = [ovTh_start ovTh_end];
    
    %% ARTEFACTS UNDER THE THRESHOLD
    unTh       = or(chanLFP<thL, [0 chanGra<thL]);
    unTh_diff  = diff(unTh);
    unTh_start = find(unTh_diff ==  1)';
    unTh_end   = find(unTh_diff == -1)';
    
%     unTh       = chanGra<thL;
%     unTh_diff  = diff(unTh);
%     unTh_start = [unTh_start; find(unTh_diff ==  1)'];
%     unTh_end   = [unTh_end  ; find(unTh_diff == -1)'];
    
    if length(unTh_end) == length(unTh_start)-1 % artefact extends the snippet
        unTh_end = [unTh_end; lenLFP];
    end
    
    if length(unTh_end)-1 == length(unTh_start) % artefact precedes the snippet
        unTh_start = [1; unTh_start];
    end
    
    unThArt    = [unTh_start unTh_end];
    
    %% START AND END OF ALL ARTEFACTS (UNDER AND OVER THE TRESHOLD)
    allArt     = [ovThArt; unThArt];                % START AND END
    artSnip    = [allArt(:,1)-75 allArt(:,2)+75];   % TAKE OUT 50ms BEFORE START AND AFTER END OF ARTEFACT
    artSnip(artSnip>lenLFP) = lenLFP;               % AFTER WINGS ARTEFACT CANNOT EXCEED  LFP SNIPPET
    artSnip(artSnip<1) = 1;                         % AFTER WINGS ARTEFACT CANNOT PRECEED LFP SNIPPET
    
    %% LOGICAL INDEX OF ARTEFACTS
    artIdx = zeros(1, lenLFP);
    for art = 1 : size(artSnip,1)
        artIdx(artSnip(art,1):artSnip(art,2)) = 1;
    end
    
    %% OPTIONAL PLOTTING
    if mkPlot == 1
        
        %% VISUALISATION
        figure(1); clf;
        hold on;
        plot(chanLFP, 'linew', 1.5, 'color', 'k');         % LFP
        plot([0 length(chanLFP)], [thH thH], 'color', 'r') % UPPER THRESHOLD
        plot([0 length(chanLFP)], [thL thL], 'color', 'r') % LOWER THRESHOLD
        
        %% MARK ARTEFACTS IN RED
        plotYlim = get(gca,'Ylim');
        for art = 1 : size(artSnip,1)
            rectangle('Position', [artSnip(art,1) plotYlim(1) artSnip(art,2)-artSnip(art,1) plotYlim(2)-plotYlim(1)], 'FaceColor', [1 0 0 0.1], 'EdgeColor', [1 1 1])
        end
        plot(artIdx, 'linew', 2, 'color', 'b') % LOGICAL INDEX OF ARTEFACTS
        
        %         subplot(212); hold on;
        %         plot(chanGra,  'k--', 'linew', 1.5);
        %         plot([0 length(chanLFP)], [thH thH], 'color', 'r')
        %         plot([0 length(chanLFP)], [thL thL], 'color', 'r')
        pause(1)
    end
    
    %% DOES LFP EXCEED TH?
    isAR(chan,:) = artIdx;
    
end
isAR = repmat(isAR, 1, 1, 100);
isAR = permute(isAR, [1 3 2]);
isAR = logical(isAR);
end % END OF FUNCTION