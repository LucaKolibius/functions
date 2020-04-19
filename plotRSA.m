function plotRSA(RSAhits, RSAmiss, comp, hitsIdx, missIdx)
% plotting
% hits
figure(666)
subplot(1,2,1)
imagesc(RSAhits)
title(['Hits ', comp]);
xlabel('Encoding');
xticks([10:10:size(hitsIdx,1)]);
ylabel('Retrieval');
set(gca,'xaxisLocation','top')
pbaspect([1 1 1]);
colorbar

% miss
if RSAmiss ~= 0
subplot(1,2,2)
imagesc(RSAmiss)
title(['Miss - ', comp]);
xlabel('Encoding');
xticks([2:2:size(missIdx,1)]);
ylabel('Retrieval');
set(gca,'xaxisLocation','top')
pbaspect([1 1 1]);
colorbar
else 
    subplot(1,2,2)
    title('No misses!');
end

set(gcf,'Position',get(0,'Screensize'));