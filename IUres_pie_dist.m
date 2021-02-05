colNU = [0.1059, 0.6196, 0.4667]; % green - NU
colIU = [0.4588, 0.4392, 0.7020]; % purple - IU
colDU = [0.8510, 0.3725, 0.0078]; % double IU

%% GU + SU
figure(1); clf;
subplot(1,3,1);
pieHand = pie([162, 5, 458], [true, true, false], {'Index Units (162)', 'Index Units (cons.) (5)', 'Single Units (458)'});
ax = gca;
ax.Colormap = [colIU; colDU; colNU];
set(findobj(pieHand, '-property', 'FaceAlpha'), 'FaceAlpha', 0.7);

% TEXT IU
pieHand(2).FontSize = 20;
pieHand(2).FontWeight = 'bold';
% pieHand(2).Position = [-0.8781 0.8179 0];

% TEXT DU
pieHand(4).FontSize = 20;
pieHand(4).FontWeight = 'bold';
pieHand(4).Position = [-0.9905 -0.2461 0];

% TEXT SU
pieHand(6).FontSize = 20;
pieHand(6).FontWeight = 'bold';
pieHand(6).Position = [0.4197 -0.9811 0];


subplot(2,3,2:3); hold on;
histHand = histogram(permGU, 70:2:150, 'Normalization', 'probability');
plot([162, 162], get(gca,'Ylim'), 'color', [0 0 0], 'linew', 5)
histHand.FaceColor = colIU;
ylabel('Probability')
yticks('')
xlabel('Number of Index Units under H_0')
xticks([70:20:170])
mAx = gca;
mAx.FontWeight = 'bold';
mAx.FontSize   = 20;
mAx.YColor = 'none';
title('p < 0.01 ')

subplot(2,3,5:6); hold on;
histHand = histogram(permIU, 0:1:10, 'Normalization', 'probability');
plot([5, 5], get(gca,'Ylim'), 'color', [0 0 0], 'linew', 5)
plot([5, 5], get(gca,'Ylim'), 'color', [0 0 0], 'linew', 5)
histHand.FaceColor = colDU;
ylabel('Probability')
yticks('')
xlabel('Number of Index Units under H_0')
xticks([0:2:8])
mAx = gca;
mAx.FontWeight = 'bold';
mAx.FontSize   = 20;
mAx.YColor = 'none';
title('p = 0.01')

