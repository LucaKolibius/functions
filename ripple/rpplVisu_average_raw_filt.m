figure(1); clf;
subplot(211)
lfpDat = rppl.lfp;
plot(mean(lfpDat,1), 'linew', 2, 'color', [0 0 0])
xlim([0 1000])
ylim([-7.5 6])
xticks('')
yticks([-4 0 4])
ylabel('\muV');
title('Unfiltered EEG (average)')
mAx = gca;
mAx.FontWeight = 'bold';
mAx.FontSize = 20;

subplot(212)
filtDat = rppl.filt;
plot(mean(filtDat,1), 'linew',2, 'color', [0 0 0])
xlabel('Time [ms]')
xlim([0 1000])
ylim([-4 4])
xticks([0:100:1000])
xticklabels([-500:100:500])
ylabel('\muV');
title('EEG from 80-140hz (average)')
mAx = gca;
mAx.FontWeight = 'bold';
mAx.FontSize = 20;


