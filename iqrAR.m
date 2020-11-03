function [isAR] = iqrAR(LFP)

numChan = size(LFP,1);    % number of channel
isAR = zeros(numChan,1);  % preallocate output
gradLFP = diff(LFP,1,2);  % difference LFP

for chan = 1 : numChan
    chanLFP = LFP(chan,:);
    chanGra = gradLFP(chan,:);
    thH = median(chanLFP) + 3*iqr(chanLFP); % iqr = prctile(x,75) - prctile(x,25)
    thL = median(chanLFP) - 3*iqr(chanLFP); % iqr = prctile(x,75) - prctile(x,25)
    
    %% DOES LFP EXCEED TH?
    if any(or (chanLFP>thH, chanLFP<thL))
        isAR(chan,1) = 1;
    end
    
    %% DOES GRADIENT EXCEED TH?
    if any(or (chanGra>thH, chanGra<thL))
        isAR(chan,1) = 1;
    end
end

end % END OF FUNCTION