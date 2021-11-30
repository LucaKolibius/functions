% housekeeping: George mentioned that microrecordings are clipped. sub-1003
% is sampled at 1024. Maybe redo MW (including interpolation and AR) myself

% spike-ripple coincidences in IU during indexing vs. during normal trials vs. coincidences in SU/GU inside indexing bundle vs coincidences in SU/GU (outside indexing bundle)
% hits vs misses (currently only hits) in encoding and retrieval
% don't average over SU!
% differentiate between indexed trials in gray units
% ripple power / frequency / duration / timing?
% only ripples from MW instead of bundle

closeVars; clear; clc;
load('X:\Luca\ripple and IU\ripplRec.mat', 'allSU');
load('X:\Luca\data\allSbj\allTrig.mat', 'trigALL', 'trigIU');

%% compare the number of spike - ripple conincidence between index units and normal single units (excluding gray untis)
spksRppl_SU = [];
for su = 1 : size(SU,2)
    
    if isempty(SU(su).su)
        continue
    end
    
    bids     = cellfun(@(x) strcmp(x, SU(su).bidsID), {trigSU.bidsID});
    sesh     = cellfun(@(x) strcmp(x, SU(su).sesh  ), {trigSU.sesh  });
    suTrig   = SU(su).su == [trigSU.su];
    trigIdx  = and(and(bids, sesh), suTrig);

    % skip gray units (1 are gray units, 0 are non-index units)
    if trigSU(trigIdx).iu == 1
        continue
    end
    
    try
    rpplPerc                = mean(SU(su).rpplRec);                                            % the percantage of the whole recording that is covered in ripples
    spksDurRppl             = sum( SU(su).rpplRec(SU(su).spks(:)) ) / size(SU(su).spks(:), 1); % the percentage of all spikes that occur during ripples
    spksRppl_SU             = [spksRppl_SU, spksDurRppl / rpplPerc];                           % a value >1 indicates that this IU spikes more often during ripples
    catch
        continue
    end
    
end
    
spksRppl_IU = [];
for su = 1 : size(IU,2)      
    rpplPerc                = mean(IU(su).rpplRec);                                            % the percantage of the whole recording that is covered in ripples
    spksDurRppl             = sum(IU(su).rpplRec(IU(su).spks)) / size(IU(su).spks, 1);         % the percentage of all spikes that occur during ripples
    spksRppl_IU             = [spksRppl_IU, spksDurRppl / rpplPerc];                           % a value >1 indicates that this IU spikes more often during ripples
end

%% ripples in indexed vs non-indexed trials vs. non-indexing bundles
% ripples in indexed trials
preCueOnly = 1;
for su = 1 : size(IU,2)
    
    % ONLY COUNT EVERY BUNDLE ONCE
        if su > 1
            if      and(strcmp( trigIU(su-1).bidsID,            trigIU(su).bidsID            ), ...
                    and(strcmp( trigIU(su-1).sesh,              trigIU(su).sesh              ), ...
                        strcmp( trigIU(su-1).wirename(1:end-1), trigIU(su).wirename(1:end-1) )))
                
                    continue
            
            end
        end
    
    % DEFINE TRIAL TIMEPOINTS
    encTrigger = trigIU(su).encTrigger(trigIU(su).hitsIdx,[1 3]);
    encTrigger = round(encTrigger * 1000);

    % TAKE ONLY THE PERIOD 1s BEFORE THE CUE
    if preCueOnly == 1
        encTrigger = [encTrigger(:,1)-1000 encTrigger(:,1)]; % precue period
    end
    
    % trial length
    trlLen     = (encTrigger(:,2) - encTrigger(:,1)) / 1000; % length of each trial in seconds
    trlLenIdx  = trlLen( logical(trigIU(su).idxTrl),:);
    trlLenNdx  = trlLen(~logical(trigIU(su).idxTrl),:);

    % WHICH TRIALS WERE INDEXED
    encTriggerIdx = encTrigger( logical(trigIU(su).idxTrl),:);
    encTriggerNdx = encTrigger(~logical(trigIU(su).idxTrl),:);
    
    % LOOP THROUGH INDEXED TRIALS
    for ind = 1:size(encTriggerIdx,1) % count how many ripples occur during each indexed trial  
        % number of ripples
        countRppls(ind) = sum(IU(su).rpplRec( encTriggerIdx(ind,1) : encTriggerIdx(ind,2) )) / trlLenIdx(ind);
        
        % spike - ripple coincidence
        rpplRec         = IU(su).rpplRec( encTriggerIdx(ind,1) : encTriggerIdx(ind,2));
        spks            = IU(su).spks   ( IU(su).spks>= encTriggerIdx(ind,1) & IU(su).spks <= encTriggerIdx(ind,2)) - encTriggerIdx(ind,1);
        
        rpplProp        = mean(rpplRec);                    % percentage of recording that has ripples
        coinc           = sum(rpplRec(spks)) / size(spks,1); % percentage of spikes that occur during ripple
        coinc_norm(ind) = coinc / rpplProp;                  % normalized coincidences
    end
    rpplNum (su,1) = nanmean(countRppls);
    coincALL(su,1) = nanmean(coinc_norm);
    
    % LOOP THROUGH NON-INDEXED TRIALS
    clear countRppls coinc_norm
    for ind = 1:size(encTriggerNdx,1) % count how many ripples occur during each non-indexed trial
        countRppls(ind) = sum(IU(su).rpplRec( encTriggerNdx(ind,1) : encTriggerNdx(ind,2) )) / trlLenNdx(ind);
        
        % spike - ripple coincidence
        rpplRec         = IU(su).rpplRec( encTriggerNdx(ind,1) : encTriggerNdx(ind,2));
        spks            = IU(su).spks   ( IU(su).spks>= encTriggerNdx(ind,1) & IU(su).spks <= encTriggerNdx(ind,2)) - encTriggerNdx(ind,1) + 1;
        
        rpplProp        = mean(rpplRec);                     % percentage of recording that has ripples
        coinc           = sum(rpplRec(spks)) / size(spks,1); % percentage of spikes that occur during ripple
        coinc_norm(ind) = coinc / rpplProp;                  % normalized coincidences
    end
    
    rpplNum (su,2) = nanmean(countRppls);
    coincALL(su,2) = nanmean(coinc_norm);

end

%% RIPPLES IN NON-INDEXING BUNDLES
rpplNum_NoN = [];
rrplNum_Gry = [];

coincSU = [];
coincGU = [];

for su = 1 : size(SU,2)
    clear countRppls
    
    % SKIP IF WE DON'T HAVE THE LFP
    if isempty(SU(su).bidsID)
        continue
    end
    
    % SKIP BUNDLES THAT HAVE AN INDEX NEURON
    if and( contains( {trigIU.bidsID}, trigSU(su).bidsID          ), ...
       and( contains( {trigIU.sesh},   trigSU(su).sesh            ), ...
            contains( {trigIU.sesh},   trigSU(su).wirename(1:end-1))))
        continue
    end
    
    % ONLY COUNT EVERY BUNDLE ONCE
        if su > 1
            if      and(strcmp( trigSU(su-1).bidsID,            trigSU(su).bidsID            ), ...
                    and(strcmp( trigSU(su-1).sesh,              trigSU(su).sesh              ), ...
                        strcmp( trigSU(su-1).wirename(1:end-1), trigSU(su).wirename(1:end-1) )))
                
                    continue
            
            end
        end
     
    % DEFINE TRIAL TIMEPOINTS
    encTrigger = trigSU(su).encTrigger(trigSU(su).hitsIdx,[1 3]);
    encTrigger = round(encTrigger * 1000);
    
    % TAKE ONLY THE PERIOD 1s BEFORE THE CUE
    if preCueOnly == 1
        encTrigger = [encTrigger(:,1)-1000 encTrigger(:,1)]; % precue period (comment out for whole trial)
    end
    
    % trial length
    trlLen     = (encTrigger(:,2) - encTrigger(:,1)) / 1000; % length of each trial in seconds
            
    % LOOP THROUGH TRIALS
    for ind = 1:size(encTrigger,1) % count how many ripples occur during each indexed trial  
        countRppls(ind) = sum(SU(su).rpplRec( encTrigger(ind,1) : encTrigger(ind,2) )) / trlLen(ind);
        
        % spike - ripple coincidence
        rpplRec         = SU(su).rpplRec( encTrigger(ind,1) : encTrigger(ind,2));
        spks            = SU(su).spks   ( SU(su).spks >= encTrigger(ind,1) & SU(su).spks <= encTrigger(ind,2)) - encTrigger(ind,1)+1;
        
        rpplProp        = mean(rpplRec);                     % percentage of recording that has ripples
        coinc           = sum(rpplRec(spks)) / size(spks,1); % percentage of spikes that occur during ripple
        coinc_norm(ind) = coinc / rpplProp;                  % normalized coincidences
    end
    
    if     trigSU(su).iu == 0
        rpplNum_NoN = [rpplNum_NoN; mean(countRppls)];
        coincSU     = [coincSU;  nanmean(coinc_norm)];
    
    elseif trigSU(su).iu == 1
        rrplNum_Gry = [rrplNum_Gry; mean(countRppls)];
        coincGU     = [coincGU;  nanmean(coinc_norm)];
    end
        
end

%% visualisation number of ripples
%  index trials
figure('units','normalized','outerposition',[0 0 1 1])
maxVal = max(rpplNum(:));
dt = 0:maxVal/10:maxVal;
subplot(221)
hold on
histogram(rpplNum(:,1), dt);
plot([mean(rpplNum(:,1)) mean(rpplNum(:,1))], get(gca, 'YLim'), 'color', 'k', 'linew', 3);
plot([median(rpplNum(:,1)) median(rpplNum(:,1))], get(gca, 'YLim'), 'color', 'k', 'linew', 3, 'linestyle', '--');
title(sprintf('Indexed Trials (M = %.2f)', mean(rpplNum(:,1))));

% non-indexed trials
subplot(223)
hold on
mhist = histogram(rpplNum(:,2), dt);
mhist.FaceColor = [0.49,0.18,0.56];
plot([mean(rpplNum(:,2))   mean(rpplNum(:,2))  ], get(gca, 'YLim'), 'color', 'k', 'linew', 3);
plot([median(rpplNum(:,2)) median(rpplNum(:,2))], get(gca, 'YLim'), 'color', 'k', 'linew', 3, 'linestyle', '--');
title(sprintf('Non-Indexed Trials (M = %.2f)', mean(rpplNum(:,2))));

% gray units 
maxVal = max(rrplNum_Gry(:));
dt = 0:maxVal/10:maxVal;
subplot(222)
hold on
mhist = histogram(rrplNum_Gry, dt);
mhist.FaceColor = [0.85,0.33,0.10];
plot([mean(rrplNum_Gry)   mean(rrplNum_Gry)  ], get(gca, 'YLim'), 'color', 'k', 'linew', 3);
plot([median(rrplNum_Gry) median(rrplNum_Gry)], get(gca, 'YLim'), 'color', 'k', 'linew', 3, 'linestyle', '--');
title(sprintf('Gray Units (M = %.2f)', mean(rrplNum_Gry)));

% normal single units
subplot(224)
hold on
mhist = histogram(rpplNum_NoN, dt);
mhist.FaceColor = [0.93,0.69,0.13];
title(sprintf('Non-Index Units (M = %.2f)', mean(rpplNum_NoN)));
plot([mean(rpplNum_NoN)   mean(rpplNum_NoN)],   get(gca, 'YLim'), 'color', 'k', 'linew', 3);
plot([median(rpplNum_NoN) median(rpplNum_NoN)], get(gca, 'YLim'), 'color', 'k', 'linew', 3, 'linestyle', '--');

% visu coincidences
% visu
dt = [1:1:25];
figure('units','normalized','outerposition',[0 0 1 1])
subplot(311)
histogram(spksRppl_IU,dt)
title(sprintf('SPK-RPPL Coincidence IU (%.2f))', mean(spksRppl_IU)))

% subplot(312)
% histogram(spksRppl_GU, dt)
% title(sprintf('SPK-RPPL Coincidence GU (%.2f))', mean(spksRppl_GU)))

subplot(313)
histogram(spksRppl_SU, dt)
title(sprintf('SPK-RPPL Coincidence SU (%.2f))', mean(spksRppl_SU)))


% spike-ripple coincidences in IU during indexing vs. during normal trials vs. coincidences in SU/GU inside indexing bundle vs coincidences in SU/GU (outside indexing bundle)  
% hits vs misses (currently only hits) in encoding and retrieval 
% don't average over SU!
% differentiate between indexed trials in gray units
% ripple power / frequency / duration / timing?
% only ripples from MW instead of bundle
