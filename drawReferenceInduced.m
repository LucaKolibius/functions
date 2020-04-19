% this belongs to xc_postClust and visualizes cluster that are either
% artefacts or reference induced (against HP wire that is referenced
% against itself)
function drawReferenceInduced(varname,x,y,z)

% "times_*" file can be named either variant
try
    load(['times_CSC_', varname(1:end-6),'.mat'], 'spikes', 'cluster_class')
catch
    load ([varname(1:end-6),'.mat'], 'spikes','cluster_class');
end

clusterNum = str2double(varname(end));
indx2 = cluster_class(:,1)==clusterNum;
clusterSpikes = spikes(indx2, :);
subplot(x,y,z)
hold on
for i = 1:size(clusterSpikes,1)
    plot(clusterSpikes(i,:), 'b');
end
zentroid = mean(clusterSpikes);
plot(zentroid, 'r', 'LineWidth',2');
title(varname)
ylim([-150 150])
xlim([1  64])
xlabel(sprintf('#Spikes: %.0f', size(clusterSpikes,1)))
hold off
end