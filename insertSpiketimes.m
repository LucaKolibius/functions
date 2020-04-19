%% insertSpiketimes
% trigger: encTrigger
% spiketimescluster=spiketimescluster
% locking: cuelocked = 1; stimlocked = 2; resplocked = 3
% timeWindow: Timewindow before and after the trigger that is considered

function spiketimes_clusterTrial=insertSpiketimes(mTrigger, spikeTimesCluster, locking, timeWindow)
spiketimes_clusterTrial={};
for ix=1:size(mTrigger,1) % number of trials  
    spiketimes_clusterTrial{1,ix} = spikeTimesCluster(find(spikeTimesCluster >= mTrigger(ix, locking) + timeWindow(1) ...
        & spikeTimesCluster <= mTrigger(ix,locking) + timeWindow(2) )) - mTrigger(ix,locking);    
end
end