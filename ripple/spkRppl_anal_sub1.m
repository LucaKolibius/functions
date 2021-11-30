% This function is used to analyse ripple activity that occours pooled over
% all MW within one bundle

% you can look at the quantity of ripples (sample points with ripples normalized by recording snippet length) and coincidences 
% (how many spikes occur on ripples, look into code for details how this is
% normalized) and the exact time points of ripples (resTime)

function [resQuant, resCoinc, resTime] = spkRppl_anal_sub1 (allSU, trigALL, trigIU, tw)

bundle = 1; % just a counter

% LOOP THROUGH ALL SINGLE UNITS (IU / GU / SU)
for su = 1 : size(allSU,2)
    
    % SKIP IF WE DON'T HAVE THE LFP
    if isempty(allSU(su).bidsID)
        continue
    end
    
    % DEFINE TRIAL BORDERS
    encTrigger = trigALL(su).encTrigger(trigALL(su).hitsIdx,[1 3]);
    encTrigger = round(encTrigger * 1000);   
    
    % DIFFERENT TIME WINDOWS OF INTEREST
    if tw == 1                                                             % PRE-CUE
        encTrigger = [encTrigger(:,1)-1000 encTrigger(:,1)];
    elseif tw == 2                                                         % PERI-CUE
        encTrigger = [encTrigger(:,1)-1000 encTrigger(:,1)+1000];
    elseif tw == 3                                                         % PERI-RESP
        encTrigger = [encTrigger(:,2)-1000 encTrigger(:,2)+1000];
    end
    
    % LOG TRIAL LENGTH
    trlLen     = (encTrigger(:,2) - encTrigger(:,1)) / 1000; % length of each trial in seconds
    
    % LOG INDEXED TRIALS
    idxd       = logical(trigALL(su).idxTrl);
    
    % Repeating Bundles are a '1'
        repBundle = 0;
        if su > 1
            if      and(strcmp( trigALL(su-1).bidsID,            trigALL(su).bidsID            ), ...
                    and(strcmp( trigALL(su-1).sesh,              trigALL(su).sesh              ), ...
                    strcmp( trigALL(su-1).wirename(1:end-1),     trigALL(su).wirename(1:end-1) )))
                
                repBundle = 1;
                
            else
                repBundle = 0;
                bundle    = bundle + 0.5;
            end
        end
        
        %% GRAY UNITS!!!!
        % SU or GU that is on an IU-BUNDLE is a '1'
        IUbundle = cellfun(@(x) x(1:end-1),    {trigIU.wirename}, 'un', 0);          % names of all the bundles with an IU
        if any( and( contains( {trigIU.bidsID}, trigALL(su).bidsID           ), ...  % current su has the same bidsID as an IU
                and( contains( {trigIU.sesh},   trigALL(su).sesh             ), ...  % same session number
                and( contains( IUbundle,        trigALL(su).wirename(1:end-1)), ...  % same wirename
                trigALL(su).iu ~= 2 ))))                                             % but it is not one of the IU
            
            onIUbundle = 1;
        else
            onIUbundle = 0;
            bundle = bundle + 0.5;
        end
        
        bundle = floor (bundle);
    
    % LOOP THROUGH INDEXED TRIALS
    for trl = 1:size(encTrigger,1) % count how many ripples occur during each indexed trial
        
        % QUANTITY OF RIPPLES
        % only count each bundle once | don't count ripples on an IU bundle again for an SU bundle
        % sort into the right variable and add the identifier (1-6)
        if onIUbundle == 0 && repBundle == 0
            
            % QUANTITY OF RIPPLES
            resQuant{bundle}(trl,1) = sum(allSU(su).rpplRec( encTrigger(trl,1) : encTrigger(trl,2) )) / trlLen(trl);
            
            % TIMERESOLUTION OF RIPPLES IN TOI
            temp = allSU(su).rpplRec( encTrigger(trl,1) : encTrigger(trl,2) );
            
            if tw ~= 0
                resTime{bundle}(trl,:) = temp;
            else
                resTime = [];
            end
            
            % ADD IDENTIFIER
            if trigALL(su).iu == 2                                     % IU
                
                if idxd(trl) == 1                                      % indexed trials
                    resQuant{bundle}(trl,2) = 1;
                else                                                   % non-indexd trials
                    resQuant{bundle}(trl,2) = 2;
                end
                
            elseif trigALL(su).iu == 1                                 % GU
                if idxd(trl) == 1                                      % indexed trials
                    resQuant{bundle}(trl,2) = 3;
                else                                                   % non-indexd trials
                    resQuant{bundle}(trl,2) = 4;
                end
                
            elseif trigALL(su).iu == 0                                 % SU (don't have indexed trials)                                             % not on a bundle that also has an index unit
                resQuant{bundle}(trl,2) = 6;
            end
            
        end
        
        %% SPIKE - RIPPLE COINCIDENCES
        %  SORT SPK-RPPL COINCIDENCE
        
        rppl                 = allSU(su).rpplRec( encTrigger(trl,1) : encTrigger(trl,2));
        spks                 = allSU(su).spks   ( allSU(su).spks>= encTrigger(trl,1) & allSU(su).spks <= encTrigger(trl,2)) - encTrigger(trl,1)+1;
        
        rpplProp             = mean(rppl);                     % percentage of recording that has ripples
        coinc                = sum(rppl(spks)) / size(spks,1); % percentage of spikes that occur during ripple
        coinc_norm           = coinc / rpplProp;               % normalized coincidences
        resCoinc{su}(trl,1)  = coinc_norm;
        
        if trigALL(su).iu == 2                                     % IU
            
            if idxd(trl) == 1                                      % indexed trials
                resCoinc{su}(trl,2) = 1;
            else                                                   % non-indexd trials
                resCoinc{su}(trl,2) = 2;
            end
            
        elseif trigALL(su).iu == 1                                 % GU
            if idxd(trl) == 1                                      % indexed trials
                resCoinc{su}(trl,2) = 3;
            else                                                   % non-indexd trials
                resCoinc{su}(trl,2) = 4;
            end
            
        elseif trigALL(su).iu == 0                                 % SU (don't have indexed trials)
            if onIUbundle == 1                                     % on a bundle that has an index unit
                resCoinc{su}(trl,2) = 5;
            else                                                   % not on a bundle that also has an index unit
                resCoinc{su}(trl,2) = 6;
            end
            
        end
               
    end % END OF TRIAL LOOP
end % END OF SU LOOP
end % END OF FUNCTION



