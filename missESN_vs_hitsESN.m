% If we were to draw a number of random trials from each single neuron. Let
% the number of draws be equal to the number of misses in that session. How
% many ESNs would we draw? -> significantly fewer than what we empirically
% find (51)

clear

load('\\analyse4.psy.gla.ac.uk\project0309\Luca\ESN_code_data\code_ESN\data\allSpks.mat', 'allSpks'); % 51 miss ESNs

nperm = 10000;
numESN = zeros(nperm,1);
for perm = 1:nperm
    for su = 1:length(allSpks)
        numMiss   = sum(allSpks(su).himiDx == 3); % number of real misses (0/2)
        reinstTrl = allSpks(su).reinstTrl(allSpks(su).himiDx == 1);        % index of which hits are reinstated trials
        
        % this is not optimal. sometimes we have more misses than hits.
        %         if length(reinstTrl)<numMiss
        %             disp(su)
        %             numMiss = length(reinstTrl);
        %         end
        %
        %         idx = logical(randi(2, [numMiss, 1])-1); % logical index
        
        % bootstrap
        mDraw = [];
        for boot = 1:numMiss
            mDraw(boot) = reinstTrl(randi(length(reinstTrl), [1 1]));
        end
            
%         if any(reinstTrl(idx))
%             numESN(perm) = numESN(perm) + 1;
%         end
        
        if any(mDraw)
            numESN(perm) = numESN(perm) + 1;
        end
        
    end % END OF SU
end % END OF PERM

figure(1); clf; hold on;
histogram(numESN, 'normalization', 'probability')
title(num2str(prctile(numESN,5)));
plot([12 12], get(gca, 'Ylim'), 'color', 'r', 'linew', 3)
% plot([36 36], get(gca, 'Ylim'), 'color', 'r', 'linew', 3)
title([num2str(mean(numESN <= 12)), ' probability that we have more miss-ESN (1/2) compared to hit-ESNs'])
% title([num2str(mean(numESN <= 36)), ' probability that we have more half-hit-ESN (1/2) compared to hit-ESNs'])
% EMPIRICALLY WE FIND 48 MISS ESNs (12 full miss, 36 halfHit)
prctile(numESN,95)

%% TemT
clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\favLibrary\data\allSpks.mat') % 27 miss ESNs

nperm = 10000;
numESN = zeros(nperm,1);
for perm = 1:nperm
    disp(perm)
    for su = 1:length(allSpks)
        numMiss   = size(allSpks(su).encTrigger,1) - sum(allSpks(su).hitsIdx,1); % number of misses
        reinstTrl = allSpks(su).reinstTrl;                                        % index of which hits are reinstated trials
        
        % this is not optimal. sometimes we have more misses than hits.
        %         if length(reinstTrl)<numMiss
        %             disp(su)
        %             numMiss = length(reinstTrl);
        %         end
        %
        %         idx = logical(randi(2, [numMiss, 1])-1); % logical index
        
        % bootstrap
        mDraw = [];
        for boot = 1:numMiss
            mDraw(boot) = reinstTrl(randi(length(reinstTrl), [1 1]));
        end
            
%         if any(reinstTrl(idx))
%             numESN(perm) = numESN(perm) + 1;
%         end
        
        if any(mDraw)
            numESN(perm) = numESN(perm) + 1;
        end
        
    end % END OF SU
end % END OF PERM

figure(1); clf; hold on;
histogram(numESN, 'normalization', 'probability')
title(num2str(prctile(numESN,5)));
plot([27 27], get(gca, 'Ylim'), 'color', 'r', 'linew', 3)

% EMPIRICALLY WE FIND 51 MISS ESNs, WHICH IS MORE THAN WE WOULD EXPECT BY
% CHANCE (SIGNIFICANTLY MORE THAN HITS)