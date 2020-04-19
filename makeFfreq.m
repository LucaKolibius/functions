%% makes a firing frequency plot under the Rasterplot
% n: histogram from makeRasterplot
% dt: binsize
% timeWindow: Timewindow before and after the trigger that is considered
% position: 17 for encCueLcoked

function makeFfreq(n,dt,timeWindow,position)
subplot(13,3,position, 'align');
fr = sum(n,1)./size(n,1)./0.175; % I have a bin each 175ms
plot(dt,fr,'ks-','LineWidth',3, 'Color', 'r'); %
%plot([0 0],[min(fr) max(fr)],'r','LineWidth',3); % red line at t=0
axis tight;
ylabel('Spikes / s', 'Color', 'k');
xlim([timeWindow(1)+0.5 timeWindow(2)-0.5]);
xlabel('Time [ms]');
set(gca, 'YLim', [0, max(fr(1:size(fr,2)-1))*1.2]);
end
