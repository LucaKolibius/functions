%% this function finds significant correlations in the hit within trial (WT) RSA for CueCue in the hippocampus
% for between trial look at fndSignBT
% for all SU time series between encoding and retrieval look at encRetTSall
clear
fetch_allSubj
subjID = allSubj{2};
try
    cd X:/Luca/data
catch
    cd /media/ldk898/rds-share/Luca/data
end

mSubject = subjID(1:end-3);
mSession = subjID(end-1:end);

% if the session name is called 1b then this line prevents an error during cd
mSubject(regexp(mSubject,'_')) = [];
if isempty(regexp(mSession,'S', 'ONCE'))
    mSession = ['S', mSession];
end

cd(mSubject)
cd(mSession)
abc = dir;
cd(abc(3).name)
cd advancedAnalysis\RSA\oldSize\hippocampus
abc = dir('RSAmat_hipp*');
load(abc(1).name, 'RSAhits_CueCue_h')

% load tables for enc and ret
cd ../../../elecLoc/spikenumbers
abc = dir('spikeNumber_hipp_*');
mCD = cd;
load(abc.name, 'encSpikeNumber_cueLocked_h', 'retSpikeNumber_cueLocked_h')
absENC = cell2mat(encSpikeNumber_cueLocked_h{3:end, :}); % needed?
absRET = cell2mat(retSpikeNumber_cueLocked_h{3:end, :}); % needed?

%% get RSA data
maskWT = logical(eye(size(RSAhits_CueCue_h,1))); % create a mask for within  trial hits
maskBT  = maskWT==0;                             % create a mask for between trial hits

WThit = RSAhits_CueCue_h(maskWT);                % extract data for within  trial from RSA
BThit = RSAhits_CueCue_h(maskBT);                % extract data for between trial from RSA
pthresh = 0.05;

%% calculate threshold
% BT distribution
% one sided z score for 0.05 = 1.645
BThit_mean = mean(BThit);
BThit_std = std(BThit);
thresh = BThit_mean+BThit_std*1.645;
idxSignWT = find(WThit>thresh); % these trial numbers are abova average similarity

%% standardize absENC + absRET to normENC + normRET
% temean = mean([absENC, absRET]  ,2);
% temstd = std ([absENC, absRET],0,2);
%
% normENC = [];
% normRET = [];
% for SU = 1 : size(absENC,1)
%     SUenc    = absENC(SU,:);
%     SUret    = absRET(SU,:);
%
%     normSUenc = (SUenc - temean(SU)) / temstd(SU);
%     normSUret = (SUret - temean(SU)) / temstd(SU);
%
%     normENC = [normENC; normSUenc];
%     normRET = [normRET; normSUret];
% end
%
% % hit miss
% himi = encSpikeNumber_cueLocked_h{2,:};
% hitIdx = strfind(himi,'hit');
%
% hitnum = 0;
% for trl = 1:size(normENC,2)
%     if isempty(hitIdx{trl}) % not a hit
%         continue % only hits
%     end
%     hitnum = hitnum+1;
%     ar1 = normENC(:,trl);
%     ar2 = normRET(:,trl);
%     temp = corrcoef(ar1,ar2); %
%     mycor(hitnum,1) = temp(2);
% end
% scatter(ar1,ar2)
allTrials = size(encSpikeNumber_cueLocked_h,2);
[tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_h, retSpikeNumber_cueLocked_h, 'hit', 0); % extracts the absolute spikenumbers for encoding and retrieval. CAREFUL: only hits and ordered by ff/pp/fp
[~, normENC, normRET] = normSpikeNumber(tempCell_enc_h, tempCell_ret_h); % normalizes spikenumber and correlates enc+ret for RSA

for indx = 1:size(idxSignWT,1) % loop over all trials in which we think there is an index
    
    % a = cell2mat(tempCell_enc_h(:,10));
    % b = cell2mat(tempCell_ret_h(:,10));
    % c = corrcoef(a,b) % 0.9469 prestandardized
    % a(end) = []; % its a zero
    % b(end) = []; % its a zero
    % c = corrcoef(a,b) % 0.9456 %% but no difference if i kick it out of the standardized matrix
    
    
    % check out which SU makes the correlation deterioate most (abs or rel doesnt matter)
    figure(idxSignWT(indx))
    
    % original correlation with all SU
    mya = normENC(:,idxSignWT(indx));
    myb = normRET(:,idxSignWT(indx));
    temp = corrcoef(mya,myb);
    refCor = temp(2);
    
    diffCor = zeros(size(mya,1),2);
    for ike = 1:size(mya,1)
        newA = mya;
        newB = myb;
        
        newA(ike) = [];
        newB(ike) = [];
        
        temp = corrcoef(newA,newB);
        diffCor(ike,1) = ike;
        diffCor(ike,2) = refCor - temp(2);
        %     newCor (ike,1) = temp(2);
    end
    
    highDec = sortrows(diffCor, 2, 'descend');
    subplot(2,5,8)
    plot(highDec(:,2), 'linew', 2) % "scree-plot" of decriment (sorted) after each SU is taken away
    axis tight
    title('sorted decrease in correlation (individual)');
    xlabel('SU');
    axH = gca;
    axH.YAxis.FontWeight = 'bold';
    axH.XAxis.FontWeight = 'bold';
    axH.FontSize = 12;
    
    %     idx = temp(:,2) >= 0.01;
    %     signSU = sum(idx);
    % %     numFigs = ceil(signSU/5);
    %
    %     toDraw = temp(idx,1);
    %
    %     if signSU>5
    %         subDim = 5;
    %     else
    %         subDim = signSU;
    %     end
    
    % plot the 5 most important SUs
    axmin = floor(min([normRET(:, idxSignWT(indx)); normENC(:, idxSignWT(indx))]));
    axmax = ceil (max([normRET(:, idxSignWT(indx)); normENC(:, idxSignWT(indx))]));
    for op = 1:5
        subh = subplot(2,5,op);
        axis square
        hold on
        
        % plots all SU in red
        h1 = scatter(normENC(:,idxSignWT(indx)),normRET(:,idxSignWT(indx)), 'filled');
        h1.CData = [0.8500 0.3250 0.0980];
        
        normRETpru               = normRET(:, idxSignWT(indx));
        normRETpru(highDec(op,1))   = [];
        normENCpru               = normENC(:, idxSignWT(indx));
        normENCpru(highDec(op,1))   = [];
        
        % only plots the remaining date in blue
        h2 = scatter(normENCpru,normRETpru, 'filled');
        h2.CData = [0 0.4470 0.7410];
        
        %     set(gca,'ylim',get(gca,'xlim'))
        set(gca,'ylim',[axmin axmax])
        set(gca,'xlim',[axmin axmax])
        yticks = xticks;
        plot(get(gca,'xlim'), [0 0 ], 'r-', 'linew',0.5)
        plot([0 0],  get(gca,'ylim'),  'r-', 'linew',0.5)
        
        if op ==1
            title(sprintf('Oldcor %.3f, decrease (w/o #%d): %.2f ', refCor, highDec(op,1), highDec(op,2)));
        else
            title(sprintf('Decrease (w/o #%d): %.2f ', highDec(op,1), highDec(op,2)));
        end
    end
    
    subplot(2,5,[6 7])
    hold on
    bar(1, highDec(1,2),'c');
    bar(2, highDec(2,2));
    bar(3, highDec(3,2));
    bar(4, highDec(4,2));
    bar(5, highDec(5,2));
    
    % aesthetics
    title('Individual contribution to the correlation');
    xticks([1:5])
    xlabel('Single Unit');
    axH = gca;
    axH.YAxis.FontWeight = 'bold';
    axH.XAxis.FontWeight = 'bold';
    axH.FontSize = 12;
    
    %% iteratively kick out SU with highest firing rate
    
    myx = [normENC(:,idxSignWT(indx)), normRET(:,idxSignWT(indx)), (normENC(:,idxSignWT(indx)) + normRET(:,idxSignWT(indx)))];
    % myx = [myx, [1:size(mysum,1)]'];
    myx = sortrows(myx,3);
    mCorr = [refCor];
    for ig = 1:size(myx,1)-2
        myx(end,:) = []; % kick out highest correlation
        temp2 = corrcoef(myx(:,1), myx(:,2)); % calculate new correlation
        mCorr = [mCorr, temp2(2)]; % extract r
        
        %     subplot(1,2,1)
        %     hold off
        %     scatter(myx(:,1),myx(:,2),'fill')
        %     title(sprintf('NewCor: %.3f', mCorr(end)));
        %
        %     subplot(1,2,2)
        %     hold on
        %     plot(mCorr, 'r')
        %     pause(1.5)
    end
    subplot(2,5,[9 10]);
    plot(mCorr,'linew',2)
    hold on
    title('Redone correlation (highest firing out first)');
    plot(get(gca,'xlim'), [0 0 ], 'r', 'linew',0.5)
    xlim([1 size(mCorr,2)]);
    xticks([1:2:size(mCorr,2)])
    xticklabels([size(mCorr,2)+1:-2:2])
    xlabel('Remaining number of Single Units');
    
    ylabel('new correlation');
    axH = gca;
    axH.YAxis.FontWeight = 'bold';
    axH.XAxis.FontWeight = 'bold';
    axH.FontSize = 12;
    trl = idxSignWT(indx);
    set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
    saveas(gcf, [subjID, sprintf('trl%.0f',trl)], 'png');
    
    
    subjID(regexp(subjID,'_')) = '-';
    cd(mCD);
    cd ../../signWT
    for op = 1:5 % SU
        su = highDec(op,1);
        mhandle = figure;
        hold on
        
        plot(normENC(su,:), 'linew',2)
        plot(normRET(su,:), 'linew',2)
        
        hand = legend({'ENC', 'RET', ''});
        xlim([1 40])
        ylabel('Firing rate in standard deviations from the mean');
        title([subjID, sprintf('-TRL%.0f-SU%.0f', trl, su)]);
        plot([trl trl],get(gca,'ylim'), 'k--', 'linew',1)
        hand.String(3) = '';
        hold off
        
        % saving figure
        set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
        saveas(gcf, [subjID, sprintf('trl%.0f_SU%.0f',trl,su)], 'png');
        axH = gca;
        axH.YAxis.FontWeight = 'bold';
        axH.XAxis.FontWeight = 'bold';
        axH.FontSize = 14;
    end
    
end
%%

%% of trial significance
idxSignWT = find(BThit>thresh); % these trial numbers are abova average similarity

% divide number by 40
dim = 40;
for x = 1:(dim*dim)-dim
    mCol(x,1) = ceil(x/(dim-1));
    mRow(x,1) = mod(x,(dim-1));
if mRow(x,1) == 0
    mRow(x,1) = (dim-1);
end
if mRow(x,1) >= mCol(x,1)
    mRow(x,1) = mRow(x,1)+1;
end
end


% % should we even consider negative Cluster?
% corrcoef(mya,myb) % initial correlation of 0.5523
% indx1 = logical(mya<0) | logical(myb <0); % indexes SU that fire more than baseline
% mya(indx1)=[]; % kick out SU that fire less than baseline
% myb(indx1)=[];
% corrcoef(mya,myb) % recalculate correlation = > 0.983
% scatter(mya,myb) % show remaining SU in scatter plot

% look at scatterplot of trl x
trl = x;
enco = normENC(:,trl);
retr = normRET(:,trl);
figure
scatter(enco,retr, 'fill');
hold on
plot(get(gca,'xlim'), [0 0 ], 'r-', 'linew',0.5)
plot([0 0],  get(gca,'ylim'),  'r-', 'linew',0.5)
corrcoef(enco,retr)

% look at the standardized firing rate for enc and ret over trials of su 14
% with a marker at trl 18
su = 14;
trl = 18;
mhandle = figure;
hold on
plot(normENC(su,:), 'linew',2)
plot(normRET(su,:), 'linew',2)
hand = legend({'ENC', 'RET', ''});
xlim([1 40])
ylabel('Firing rate in standard deviations from the mean');
title([subjID, sprintf('-TRL%.0f-SU%.0f', trl, su)]);
plot([trl trl],get(gca,'ylim'), 'k--', 'linew',1)
hand.String(3) = '';
hold off
axH = gca;
        axH.YAxis.FontWeight = 'bold';
        axH.XAxis.FontWeight = 'bold';
        axH.FontSize = 14;