%% automatically find artefacts in the LFP based on _threshold_ std from the mean
% concatenates all trials per MW
function [ del_ts ] = LFP_artDet( data, dir_saves, threshold )

for curChan = 1:length( data.label )%loop over channels
    % concatenate trials
    lfp_data = [];
    lfp_time = [];
    for trl = 1:length(data.trial)
        lfp_data = [lfp_data, data.trial{trl}(curChan,:)];
        lfp_time = [lfp_time, data.sampleinfo(trl,1):1:data.sampleinfo(trl,2)];
    end

    % compute z-score of LFP
    m  = mean(lfp_data,2);
    sd = std(lfp_data,0,2);
    z = ( lfp_data-m )./ sd;
    z = abs(z);
    
    %%
    abvTH = z >= threshold;
    abvTH = diff(abvTH);
%     fromTo = [(find(abvTH == 1)+1)' (find(abvTH == -1))'];
    fromIdx = find(abvTH ==1)+1;
    toIdx = find(abvTH == -1);
    fromToIdx = [(find(abvTH == 1)+1)' (find(abvTH == -1))'];
    fromTo = lfp_time(fromToIdx);
    
    %%
    % keep timepoints where LFP-activity is out of range (outliers)
    del_ts{curChan} = fromTo;
end

save([dir_saves '.mat'], 'del_ts', '-v7.3');
end