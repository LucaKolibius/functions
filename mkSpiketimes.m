% create a spiketime graphic using scattered points representing each spike
% in one cluster

function h1 = mkSpiketimes(dts, it)
spikeFiles = dir('times_*.mat');
mySpikefile = spikeFiles(it).name;
load(mySpikefile, 'cluster_class');


h1 = figure(99);
dotSize = 30;
color = [0 0 0; 0 1 0; 0 0 1; 1 0 0; 0 1 1; 1 1 0; 1 0 1];

for ia = 0:max(cluster_class(:,1))
subplot(max(cluster_class(:,1))+1,1,ia+1);
hold on
title(sprintf('Cluster #%d', ia));

spiketimes = cluster_class(cluster_class==ia,2);
spiketimes = spiketimes/32000;

y = rand(size(spiketimes));

scatter(spiketimes,y,dotSize, color(ia+1,:), '.');
xlim([0 dts.end+dts.start])

plot([dts.start dts.start],[0 1],'r','LineWidth',3); % red line at exp beginning
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis

%       n = hist(x,dts.dt);
end

xlabel(['Time in Seconds', char(10),'(From Start of Experiment)']); % label x-axis
end