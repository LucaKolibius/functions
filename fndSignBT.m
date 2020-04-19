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
idxSignBT = find(BThit>thresh); % these trial numbers are above average similarity
allTrials = size(encSpikeNumber_cueLocked_h,2);
[tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_h, retSpikeNumber_cueLocked_h, 'hit', 0); % extracts the absolute spikenumbers for encoding and retrieval. CAREFUL: only hits and ordered by ff/pp/fp
[~, normENC, normRET] = normSpikeNumber(tempCell_enc_h, tempCell_ret_h); % normalizes spikenumber and correlates enc+ret for RSA

% get coordinates in the RSA of significant BT trials
dim = size(RSAhits_CueCue_h,2);
clear mCol mRow
for loop = 1:size(idxSignBT,1)
    x = idxSignBT(loop);
    mCol(loop,1) = ceil(x/(dim-1));
    mRow(loop,1) = mod(x,(dim-1));
if mRow(loop,1) == 0
    mRow(loop,1) = (dim-1);
end
if mRow(loop,1) >= mCol(loop,1)
    mRow(loop,1) = mRow(loop,1)+1;
end
end

retTrl = mRow;
encTrl = mCol;

for indx = 1:size(idxSignBT,1) % loop over all significant trials
    figure(idxSignBT(indx))
    
    % original correlation with all SU
    mya = normENC(:,encTrl(indx));
    myb = normRET(:,retTrl(indx));
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
    
    axmin = floor(min([normRET(:, retTrl(indx)); normENC(:, encTrl(indx))]));
    axmax = ceil (max([normRET(:, retTrl(indx)); normENC(:, encTrl(indx))]));
    
    for op = 1:5
        subh = subplot(2,5,op);
        axis square
        hold on
        
        % plots all SU in red
        h1 = scatter(normENC(:, encTrl(indx)),normRET(:, retTrl(indx)), 'filled');
        h1.CData = [0.8500 0.3250 0.0980];
        
        normRETpru               = normRET(:, retTrl(indx));
        normRETpru(highDec(op,1))   = [];
        normENCpru               = normENC(:, encTrl(indx));
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
    bar(1, highDec(1,2));
    bar(2, highDec(2,2));
    bar(3, highDec(3,2));
    bar(4, highDec(4,2));
    bar(5, highDec(5,2));
    
    % aesthetics
    title('Correlation decrease for each SU');
    xticks([1:5])
    xlabel('Single Unit');
    axH = gca;
    axH.YAxis.FontWeight = 'bold';
    axH.XAxis.FontWeight = 'bold';
    axH.FontSize = 12;
    
     myx = [normENC(:,encTrl(indx)), normRET(:,retTrl(indx)), (normENC(:,encTrl(indx)) + normRET(:,retTrl(indx)))];
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
    trl = idxSignBT(indx);
    set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
    saveas(gcf, ['BT_',subjID, sprintf('_tr%.0f',trl)], 'png');
    