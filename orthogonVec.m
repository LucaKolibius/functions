function [ data, success, bundles] = orthogonVec( data )

bundles = unique(cellfun(@(x)x(1:end-1),data.label, 'un',0));


% loop through labels
trials  = 1 : length(data.trial); 
success = [];

for bund = 1 : length(bundles) % LOOP THROUGH BUNDLES
    
    curMWs = contains(data.label,bundles{bund});
     
    for curTrl = trials % LOOP THROUGH TRIALS
        
        
        mwIdx = find(curMWs);
        LFP = [];
        for mw = 1:length(mwIdx) % LOOP THROUGH ALL MW ON THE CURRENT BUNDLE
            
            oneMW   = mwIdx(mw);                                % take one MW
            otherMW = mwIdx(mwIdx~=mwIdx(mw));                  % take all the other MW
            oneLFP = data.trial{curTrl}(oneMW,:);               % extract LFP from oneMW
            otherLFP = mean(data.trial{curTrl}(otherMW,:), 1);  % mean of all other chans
            
            % ORTHONORMALIZATION
            % multiply both LFPs element wise and and sum that vector up
            % multiply otherLFP with itself element wise and sum that vector up 
            % divide #1 with #2
            % use that to scale otherLFP
            
            Xp = ( sum( otherLFP.*oneLFP )./sum( otherLFP.*otherLFP ) ) * otherLFP;  % orthonormalise
            newLFP = oneLFP - Xp;                                                    % subtract
            success(mwIdx(mw),curTrl) = dot(newLFP, otherLFP);                       % check if vectors are orthogonal
            LFP(oneMW-mwIdx(1)+1,:) = newLFP;
        end
        data.trial{curTrl}(curMWs,:) = LFP;
        if(curTrl == length(data.trial))
            disp(['CHECK BUNDLE ' bundles{bund} ' VECTOR ORTHOGONALITY ' ...
                ': ' num2str(mean(mean(success(mwIdx, trials)))) ' (SHOULD BE CLOSE TO 0)'])
        end
    end
    
    
%     data.label(curMWs) = cellfun(@(x)[x ' GS REF'], data.label(curMWs), 'un', 0);
end

end

