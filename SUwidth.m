clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_encTwo_retResp_noCN.mat', 'allSpks');

for su = 1:length(allSpks)
    [~, subjID, ~, ~, ~, ~, ~, averageWS, ~, ~] = loadInDat(allSpks, su);
    sr = allSpks(su).initSR;
    
    [peakVal, peakPos]     = max(averageWS);
    [troughVal, troughPos] = min(averageWS);
    
    if abs(troughVal) > abs(peakVal)
        averageWS = averageWS*-1;
        
        [peakVal, peakPos]     = max(averageWS);
        [troughVal, troughPos] = min(averageWS);
    end
    
    
    halfMaxi = (peakVal-troughVal)/2;
    
    interpFac = 100;
    WS_upsamp = interp(averageWS, interpFac);
    
    diffHalfMaxi        = abs(WS_upsamp-halfMaxi);
    [sortDiff, sortIdx] = sort(diffHalfMaxi);
    widthStart          = sortIdx(1); % the start of the half maximum
    
    foundBoth = 0;
    count     = 2;
    while foundBoth == 0
        if abs(widthStart-sortIdx(count)) < interpFac+1
            count = count + 1;
        else
            widthEnd  = sortIdx(count);
            foundBoth = 1;
        end
    end
    
    fullWidth = abs(widthEnd-widthStart);
    fullWidth = fullWidth / (sr*interpFac); % in ms
    suClas(su).fullWidth = fullWidth;
    suClas(su).spikeHeight = halfMaxi*2;
    suClas(su).isIU = allSpks(su).iu;
    
    for period = 1:2
        if period == 1
            encRet = 'enc';
        else
            encRet = 'ret';
        end
        
        %% TRIGGER IS IN SECONDS
        switch strcmp(encRet, 'enc')
            case 0 % RETRIEVAL TRIGGER
                trig = allSpks(su).retTrigger(allSpks(su).hitsIdx, [1 3]);
            case 1 % ENCODING TRIGGER
                trig = allSpks(su).encTrigger(allSpks(su).hitsIdx, [1 3]);
        end
        
        idxTrl     = allSpks(su).idxTrl;
        
        %% SPIKES IN SECONDS
        clusterSpikes      =  allSpks(su).spks/1000;
        [spksSeg, trlLen]  =  insertSpiketimes2(trig, clusterSpikes, [1 2], [0 0]); % 3 seconds prior to cue trigger until 1 second after response trigger
        spksSeg            =  cell2mat(cellfun(@length, spksSeg, 'un', 0));
        
        fireRate    = spksSeg ./ trlLen;
        fireRateIdx = sum(spksSeg( idxTrl)) / sum(trlLen( idxTrl));
        fireRateNdx = sum(spksSeg(~idxTrl)) / sum(trlLen(~idxTrl));
        
        
        fanoFac = (std(fireRate, 0, 1)^2 ) / mean(fireRate,1);
        
        suClas(su).fullWidth = fullWidth;
        
        if period == 1 % ENCODING
            suClas(su).encFireIdx = fireRateIdx;
            suClas(su).encFireNdx = fireRateNdx;
            suClas(su).encFano    = fanoFac;
        else
            suClas(su).retFireIdx = fireRateIdx;
            suClas(su).retFireNdx = fireRateNdx;
            suClas(su).retFano    = fanoFac;
        end
    end
end

%% VISUALISATION
colNU = [0.1059, 0.6196, 0.4667]; % green - NU
colIU = [0.4588, 0.4392, 0.7020]; % purple - IU

isIU = logical([suClas.isIU]);
figure(1); clf;
subplot(211)
allWidthI = [suClas(isIU).fullWidth]; % IU
dt = 0:0.000025:0.0005;
hhist = histogram(allWidthI,dt, 'normalization', 'probability');
hhist.FaceColor = colIU;

hold on; % SU
allWidthS = [suClas(~isIU).fullWidth];
hhist = histogram(allWidthS,dt, 'normalization', 'probability');
hhist.FaceColor = colNU;

xlabel('Spike Width (FWHM) [ms]')
ylabel('Probability')
legend({'ESN','Single Units'})
legend('boxoff')
set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);
box off

[~,pWidth,~,~] = ttest2(allWidthI, allWidthS)

subplot(212)
allHeightI = [suClas(isIU).spikeHeight];
dt = 0:20:300;
hhist = histogram(allHeightI, dt, 'normalization', 'probability');
hhist.FaceColor = colIU;

hold on
allHeightS = [suClas(~isIU).spikeHeight];
dt = 0:20:300;
hhist = histogram(allHeightS, dt, 'normalization', 'probability');
hhist.FaceColor = colNU;

[~,pHeight,~,~] = ttest2(allHeightI, allHeightS)

xlabel('Spike Height [\muV]')
ylabel('Probability')
legend({'ESN','Single Units'})
legend('boxoff')
set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);
box off

%%
figure(2); clf;
subplot(211);
fanoEncI = [suClas(isIU).encFano];
dt = 0:0.2:5;
hhist = histogram(fanoEncI, dt, 'normalization', 'probability');
hhist.FaceColor = colIU;

hold on
fanoEncS = [suClas(~isIU).encFano];
dt = 0:0.2:5;
hhist = histogram(fanoEncS, dt, 'normalization', 'probability');
hhist.FaceColor = colNU;

[~,pFanoEnc,~,~] = ttest2(fanoEncI, fanoEncS)

title('Fanofactor encoding')
legend({'ESN','Single Units'})
legend('boxoff')
xlabel('Fanofactor')
ylabel('Probability')
set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);
box off

subplot(212);
fanoRetI = [suClas(isIU).retFano];
hhist = histogram(fanoRetI, dt, 'normalization', 'probability');
hhist.FaceColor = colIU;

hold on
fanoRetS = [suClas(~isIU).retFano];
dt = 0:0.2:5;
hhist = histogram(fanoRetS, dt, 'normalization', 'probability');
hhist.FaceColor = colNU;

[~,pFanoRet,~,~] = ttest2(fanoRetI, fanoRetS)

title('Fanofactor retrieval')
legend({'ESN','Single Units'})
legend('boxoff')
xlabel('Fanofactor')
ylabel('Probability')
set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);
box off

%%
figure(3); clf;
subplot(221);
fireEncReinst = [suClas(isIU).encFireIdx];
dt = 0:2.5:30;
hhist = histogram(fireEncReinst, dt, 'normalization', 'probability');
hhist.FaceColor = colIU;

xlabel('Firing rate [hz]')
ylabel('Probability')
title('Firing Rate: Encoding Reinstated')
ylim([0 1])
yticks([0:0.2:1])
set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);
box off

subplot(222);
fireEncNon = [suClas(isIU).encFireNdx];
dt = 0:2.5:30;
hhist = histogram(fireEncNon, dt, 'normalization', 'probability');
hhist.FaceColor = colIU;

hold on
fireEncSU = [suClas(~isIU).encFireNdx];
dt = 0:2.5:30;
hhist = histogram(fireEncSU, dt, 'normalization', 'probability');
hhist.FaceColor = colNU;

[~,pEnc_NonVSsu,~,~]   = ttest2(fireEncNon, fireEncSU)
[~,pEnc_RenstVSsu,~,~] = ttest2(fireEncReinst,  fireEncSU)
mean(fireEncNon)
mean(fireEncReinst)
mean(fireEncSU)


legend({'ESN','Single Units'})
set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);
legend('boxoff')
box off
xlabel('Firing rate [hz]')
ylabel('Probability')
title('Firing Rate: Encoding Non-Reinstated')
ylim([0 1])
yticks([0:0.2:1])

subplot(223);
fireRetReinst = [suClas(isIU).retFireIdx];
dt = 0:2.5:30;
hhist = histogram(fireRetReinst, dt, 'normalization', 'probability');
hhist.FaceColor = colIU;

xlabel('Firing rate [hz]')
ylabel('Probability')
title('Firing Rate: Retrieval Reinstated')
ylim([0 1])
yticks([0:0.2:1])
set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);
box off

subplot(224);
fireRetNon = [suClas(isIU).retFireNdx];
dt = 0:2.5:30;
hhist = histogram(fireRetNon, dt, 'normalization', 'probability');
hhist.FaceColor = colIU;

hold on;
fireRetSU = [suClas(~isIU).retFireNdx];
dt = 0:2.5:30;
hhist = histogram(fireRetSU, dt, 'normalization', 'probability');
hhist.FaceColor = colNU;

[~,pRet_NonVSsu,~,~]    = ttest2(fireRetNon, fireRetSU)
[~,pRet_RenstVSsu,~,~]  = ttest2(fireRetReinst, fireRetSU)
mean(fireRetNon)
mean(fireRetReinst)
mean(fireRetSU)

% % ranksum(fireIenc,fireSenc)
% % ranksum(fireIret,fireSret)
% % [~,p]=kstest(fireIenc)
% % [~,p]=kstest(fireSenc) 
% % [~,p]=kstest(fireIret) 
% % [~,p]=kstest(fireSret)


legend({'ESN','Single Units'})

xlabel('Firing rate [hz]')
ylabel('Probability')
title('Firing Rate: Retrieval Non-Reinstated')
ylim([0 1])
yticks([0:0.2:1])
set(findobj(gca,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);
legend('boxoff')
box off


figure(1);
print('\\analyse4.psy.gla.ac.uk\Project0309\Writing the Index\figures\su features\width_height','-dpng','-r300')

figure(2);
print('\\analyse4.psy.gla.ac.uk\Project0309\Writing the Index\figures\su features\fanofactor_encRet','-dpng','-r300')

figure(3);
print('\\analyse4.psy.gla.ac.uk\Project0309\Writing the Index\figures\su features\firingRate','-dpng','-r300')
