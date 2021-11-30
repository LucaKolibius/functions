% chosing optimal time windows
% this is necessary because we have nans in the lfp (previously artefacts)
% in the later script I cannot use spikes that occur in the same segment as
% nans. This is my approach to salvage more spikes!
function [tw1, tw2] = fndLFPsegs(data, idxSpks, nwin)
dat       =  data.lfp_noArt;
spikelocs =  round(idxSpks);
maxim     =  0; % number of spikes in segments not maximized yet
tw        =  3000; % initial segment length
expDur    =  size(dat,2);

while maxim == 0   
    cont = 1;
    
    tw1 = spikelocs-tw*0.5;
    tw2 = spikelocs+tw*0.5;
    
    % fix trial window if it deceeds 1 or exceeds the experiment duration
    tw2(tw1<0) = tw;
    tw1(tw1<0) = 1; 

    tw1(tw2>expDur) = expDur - tw;
    tw2(tw2>expDur) = expDur;
    
    for seg = 1:size(tw1,1)
        lostSpks = 0;
        dir = 1;
        while and(any(any(isnan(dat(:,tw1(seg) : tw2(seg))))) == 1, cont) % this is the first segment. is there a nan value in any of the channels?. takes 10ms
            % as long as there are nans anywhere, the LFP snippet will be
            % shifted left
            if and( and(spikelocs(seg)+nwin < tw2(seg), dir == 1), tw1~=1)
                tw1(seg) = tw1(seg)-1;
                tw2(seg) = tw2(seg)-1;
                maxim = 1;
                % or right
            elseif and(spikelocs(seg)-nwin >= tw1(seg) , tw2 ~= expDur)
                dir = 0;
                tw1(seg) = tw1(seg) + 1;
                tw2(seg) = tw2(seg) + 1;
                maxim = 1;
                % if it cannot find a fit, stop working with this tw
            else %if spikelocs(seg)-nwin < tw1(seg)
                cont = 0;
                maxim = 0;
            end
        end
    end
    
    % reduce trial length by 2 samples
    tw = tw-2;
    
    % minimum segment lenght is 150
    if tw == 2998 %1498
        maxim = 1;
    end
end