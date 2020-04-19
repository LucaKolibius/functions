%% insertSpiketimes
% trigger: encTrigger
% spiketimescluster=spiketimescluster
% locking: cuelocked = 1; stimlocked = 2; resplocked = 3
% timeWindow: Timewindow before and after the trigger that is considered

% in insertSpiketimes2 the function can deal with two different lockings
% for trialstart and trialfinish

function [spiketimes_clusterTrial, twlength] =insertSpiketimes2(mTrigger, spikeTimesCluster, locking, timeWindow)
spiketimes_clusterTrial={};
for ix=1:size(mTrigger,1) % number of trials  
    spiketimes_clusterTrial{1,ix} = spikeTimesCluster(find(spikeTimesCluster >= mTrigger(ix, locking(1)) + timeWindow(1) ...
        & spikeTimesCluster <= mTrigger(ix,locking(2)) + timeWindow(2) )) - mTrigger(ix,locking(1));
    
    twlength = (mTrigger(ix,locking(1)) + timeWindow(1)) - (mTrigger(ix,locking(2)) + timeWindow(2));
end
end