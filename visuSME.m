% here I visualize the results of the SME for hits and misses

mH = figure(2);
col1 = [0.552941176470588,0.827450980392157,0.780392156862745];
col2 = [0.984313725490196,0.501960784313726,0.447058823529412];

upperBound_hits = (gausConvENC_hitsALL_m + gausConvENC_hitsALL_se)';
lowerBound_hits = (gausConvENC_hitsALL_m - gausConvENC_hitsALL_se)';

upperBound_miss = (gausConvENC_missALL_m + gausConvENC_missALL_se)';
lowerBound_miss = (gausConvENC_missALL_m - gausConvENC_missALL_se)';

x = (1:size(gausConvENC_hitsALL_m, 2))';
hold on

% shade from stim to stim+1s (14625 to 19625)
colShade = [0.8275 0.8275 0.8275];
shadeWindow = 14625:1:19625;
ShadeX = [shadeWindow fliplr(shadeWindow)];
ShadeY = [zeros(1,5001)-0.99 zeros(1,5001)+5];

%% Visualization
% Shade
mShade = fill(ShadeX, ShadeY, colShade);
set(mShade,'EdgeColor',colShade); % consider getting rid of the black edges


% hitss
hP_hits = patch([x; x(end:-1:1); x(1)], [lowerBound_hits; upperBound_hits(end:-1:1); lowerBound_hits(1)], 'r');
hL_hits = line(x,gausConvENC_hitsALL_m);

set(hP_hits, 'facecolor',col1, 'edgecolor','none', 'FaceAlpha',0.4);
set(hL_hits, 'color', col1, 'marker', '.');

% misses
hP_miss = patch([x; x(end:-1:1); x(1)], [lowerBound_miss; upperBound_miss(end:-1:1); lowerBound_miss(1)], 'r');
hL_miss = line(x,gausConvENC_missALL_m);

set(hP_miss, 'facecolor',col2, 'edgecolor','none', 'FaceAlpha',0.4);
set(hL_miss, 'color',col2, 'marker','.');

% straight line at y = 0
line([0,29250],[0,0], 'LineStyle','--', 'Color','k', 'LineWidth',1.5)

% legend
legend([hP_hits hP_miss], {'Hits', 'Miss'},'FontSize',16, 'FontWeight','bold')
legend('boxoff');

% aesthetics
ylabel({'Firing Rate [\delta]';''}, 'FontSize',18, 'FontWeight','bold'); % y-axis label
yticks([-0.5 -0.25 0 0.25 0.5]);
ylim([-1 1]);
% xlabel('Time');
title({'Subsequent Memory Effect';''}, 'FontSize', 25); % title

mAx = gca;
mAx.YAxis.FontWeight = 'bold'; % y-axis bold
mAx.XAxis.FontWeight = 'bold'; % x-axis bold
mAx.XAxis.LineWidth = 2; % x-axis width
mAx.YAxis.LineWidth = 2; % y-axis width
mAx.FontSize = 18; % fontsize of xlabel
mAx.YGrid = 'on';

xticks([4625, 14625])
xlim([1 29250])
xticklabels({'Cue Onset', 'Stimulus Onset'});

% % unused
% fill([1:5001 fliplr(1:5001)], [ones(1,5001)*1, ones(1,5001)*2], 'k')
% patch([x fliplr(x)], [y1 fliplr(y2)], 'g')

