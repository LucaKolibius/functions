clear


col = [255,255,204;
161,218,180;
65,182,196;
34,94,168];
col = col/255; 

col1 = col(1,:);
col2 = col(2,:);
col3 = col(3,:);
col4 = col(4,:);

% colNU = [0.6980 0.8745 0.5412]; % green - NU
% colIU = [0.1216 0.4706 0.7059]; % purple - IU
% colDU = [0.6510 0.8078 0.8902]; % double IU



%% GU + SU
figure(1); clf;
% subplot(1,3,1);
% pieHand = pie([162, 5, 458], [true, true, false], {'Episode Specific Neuron (162)', 'Episode Specific Neuron (two cat.) (5)', 'Single Units (458)'});
% pieHand = pie([resGU.num-resIU.num, resIU.num, size(allSpks,2)-resGU.num], [true, true, true], {'', '', ''});
pieHand = pie([478 147], [true, true], {'', ''});
ax = gca;
ax.Colormap = [col2; col4];
set(findobj(pieHand, '-property', 'FaceAlpha'), 'FaceAlpha', 1);

% % legend({'Single units', 'ESN (without CC)', 'ESN (FF / PP)', 'ESN'}, 'FontSize',25, 'FontWeight', 'normal');
% legend('boxoff');


% % TEXT IU
% pieHand(2).FontSize = 20;
% pieHand(2).FontWeight = 'bold';
% % pieHand(2).Position = [-0.8781 0.8179 0];
% 
% % TEXT DU
% pieHand(4).FontSize = 20;
% pieHand(4).FontWeight = 'bold';
% pieHand(4).Position = [-0.9905 -0.2461 0];
% 
% % TEXT SU
% pieHand(6).FontSized = 20;
% pieHand(6).FontWeight = 'bold';
% pieHand(6).Position = [0.4197 -0.9811 0];

%%
% subplot(2,3,2:3); hold on;
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp.mat', 'allSpks','permIU', 'permGU', 'resIU', 'resGU');
figure(2); clf; %hold on;
subplot(131); hold on;
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp.mat')
histHand = histogram(permGU, min(permGU):2:max(permGU), 'Normalization', 'probability');
ylim([0 0.1])
yticks([0 0.04 0.08])
%         xlim([70 200])
% ylabel('Probability')

plot([resGU.num, resGU.num], get(gca,'Ylim'), 'color', [0 0 0], 'linew', 3)
histHand.FaceColor = col1;
histHand.FaceAlpha = 1;
% yticks('')
% yticklabels({'0', '0.02', '0.04', '0.06', '0.08', '0.1'})
xlabel('Number of ESN under H_0')

xticks([70:20:170])
% mAx.YColor = 'none';
% title('p = .0145 ')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(gca, 'YAxisLocation', 'right');
view([90 90])
xlim([70 170])

%%
% subplot(2,3,5:6); hold on;
% figure(3); clf; hold on;
subplot(133); hold on;
ylim([0 0.35])
histHand = histogram(permIU, min(permIU):1:max(permIU), 'Normalization', 'probability');
plot([resIU.num, resIU.num], get(gca,'Ylim'), 'color', [0 0 0], 'linew', 3)
histHand.FaceColor = col3;
histHand.FaceAlpha = 1;

% ylabel('Probability')
yticks([0 0.15 0.3])
xticks([0:2:8])
xlim([-2.5 9.5])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(gca, 'YAxisLocation', 'right');
view([90 90])
% mAx.YColor = 'none';
% title('p = .0324')

%%
subplot(132); hold on;
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp_noCN.mat')
histHand = histogram(permGU, min(permGU):2:max(permGU), 'Normalization', 'probability');
plot([resGU.num, resGU.num], get(gca,'Ylim'), 'color', [0 0 0], 'linew', 3)
histHand.FaceColor = col4;
histHand.FaceAlpha = 1;

ylabel('Probability')
yticks([0 0.05 0.1 0.2 0.3])
% xlabel('Number of Index Units under H_0')
% xlim([-0.5 9.5])
xticks([40:20:110])
xlim([40 150])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
set(gca, 'YAxisLocation', 'right');
view([90 90])

%% NEW NOT TILTED
%%
figure(2); clf ; hold on;
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp_noCN.mat')
histHand = histogram(permGU, min(permGU):2:max(permGU), 'Normalization', 'probability');
histHand.FaceColor = col4;
histHand.FaceAlpha = 1;

ylabel('Probability')
xlabel('Number of ESN under H_0')

yticks([0 0.05 0.1])
xticks([80:20:140])
xlim([70 150])
ylim([0 0.1])
plot([resGU.num, resGU.num], get(gca,'Ylim'), 'color', [0 0 0], 'linew', 3)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15,'FontWeight','Bold', 'LineWidth', 2);
