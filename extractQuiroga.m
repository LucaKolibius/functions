% nun = cellfun(@(x) mean(x,1), allNunSU, 'un', 0);
% nun = vertcat(nun{:});
% nunM = nanmean(nun,1);
% nunS = nanstd(nun,0,1) / sqrt(size(nun,1));
%
%
% tun = cellfun(@(x) mean(x,1), allTunSU, 'un', 0);
% tun = vertcat(tun{:});
% tunM = nanmean(tun,1);
% tunS = nanstd(tun,0,1) / sqrt(size(tun,1));
%
% figure(1); clf; hold on;
%
% plotDT  = [0:1:1000];
% shadeDT = [plotDT, fliplr(plotDT)];
%
% plot(plotDT, nunM, 'color', [0 0 0], 'linew', 2.5);
% inBetweenNun = [nunM-nunS, fliplr(nunM+nunS)];
% fillHandNun  = fill(shadeDT, inBetweenNun, [0 0 0]);
% fillHandNun.FaceAlpha = 0.25;
%
% plot(plotDT, tunM, 'color', [1 0 0], 'linew', 2);
% inBetweenTun = [tunM-tunS, fliplr(tunM+tunS)];
% fillHandTun  = fill(shadeDT, inBetweenTun, [1 0 0]);
% fillHandTun.FaceAlpha = 0.25;
%
%
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',35,'FontWeight','Bold', 'LineWidth', 2);
% xlabel('Time [ms]')
% ylabel('Standardized Firing Rate')
% xlim([0 1000])
% ylim([-2 6])
%
% for tp = 1:50:size(tun,2)
% scathandler = scatter(ones(size(tun,1),1)*tp, tun(:,tp), 6, [0 0 0]);
% scathandler.MarkerEdgeAlpha = 0.5;
% scathandler.MarkerFaceAlpha = 0.5;
% end
%
% L(1) = plot(nan, nan, 'color', [0 0 0]);
% L(2) = plot(nan, nan, 'color', [1 0 0]);
% legend(L, {'Non-tuned trials', 'Tuned-trials'}, 'FontSize',25, 'FontWeight', 'bold');
% legend('boxoff')




%% extracts spiketimes from CN in declans data
function extractQuiroga
cd('\\analyse4.psy.gla.ac.uk\project0309\Declans_Tuning_Data\tuning');
allSbj = dir('sub-00*');

allTun = [];
skipped  = [];

% allTunSU    = [];
% allNunSU    = [];

dt      = -0.625:0.001:1.125+0.001;
dt      = dt-0.0005; % center around 0

mlength = 251;
mKernel = gausswin(mlength);

%% gaussian kernel
% mlength = [-0.075:0.0002:0.075];
% mSigma  = 0.02; % 20ms
% mKernel = normpdf(mlength,0,mSigma);
% mKernel = mKernel/max(mKernel); % normalize peak to 1
% mlength = 251;
% mKernel = gausswin(mlength);
%
% dt      = -0.625:0.001:1.125+0.001;
% dt      = dt-0.0005; % center around 0


for sbj = 1:size(allSbj,1)
    cd([allSbj(sbj).folder, filesep, allSbj(sbj).name]);
    
    allSesh = dir('20*');
    
    for sesh = 1:size(allSesh,1)
        %         allSU = {};
        
        if strcmp(allSesh(sesh).name(end), 'U') % U session's are unusable
            continue
        end
        
        try
            cd([allSesh(sesh).folder, filesep, allSesh(sesh).name]);
            cd('SU Analysis')
            load('allSUdata', 'allCl', 'ttl_vec', 'image_cell');
            load('allRespData', 'resp_mat');
        catch
            skipped = [skipped; sbj, sesh];
            continue
        end
        numSU    = size(allCl,1);
        
        for su = 1:numSU
            allSU = insertSpiketimes(ttl_vec, allCl{su,2}, 1, [0.300 1.00]); % 700ms
            allBL = insertSpiketimes(ttl_vec, allCl{su,2}, 1, [-0.500 -0.200]); % 300ms
            allFir = insertSpiketimes(ttl_vec, allCl{su,2}, 1, [-0.625 1.125]); % for the firing rate analysis
            
            
            allSU = cellfun(@length, allSU, 'un', 0); allSU = cell2mat(allSU); % trial period
            allBL = cellfun(@length, allBL, 'un', 0); allBL = cell2mat(allBL); allBL = allBL./300.*700; % baseline
            
            % check if SU is tuned to any image
            for img = 1:size(image_cell,1)
                idx = image_cell{img,2};
                curMed = median(allSU(idx));
                blMean = mean(allBL);    % mean of baseline / all pictures
                blStd  = std(allBL,0,2); % std  of baseline / all pictures
                firTh  = blMean + 1*blStd; % SET BACK TO 5!!!
                
                tunImg = false(size(image_cell,1),1);
                if curMed >= 2 && curMed > firTh
                    image_cell{img,2+su} = 1;
                    tunImg(img) = 1;
                else
                    image_cell{img,2+su} = 0;
                end
            end
                        
            tunTrials = [image_cell{tunImg,2}];
            
            if sum(tunImg) == 0 % no tuning to any image
                continue
            end
            
            % CONVOLVE ALL TRIALS TO SPIKE SERIES
            allConvRaw = [];
            for trl = 1:size(allFir,2)
                [x,~]      = histcounts(allFir{trl}, dt);
                spkConv    = conv(x, mKernel, 'same');
                allConvRaw = [allConvRaw; spkConv];
                
            end
            allConv = allConvRaw(:,126:end-125); % cut off wings
            
            BL = mean(allConv(:,1:500),2); % baseline is 500ms
            allConv = allConv(:,501:end);  % trial data
            %             BL = mean(allConv,2); % over time // whole trial! (old)
            BLmean = mean(BL); % mean for z-scoring
            BLstd  = std(BL, 0); % std for z-scoring
            BLstd  = BLstd + 0.1; % add 0.1 for sparse SU
            
            % z-scoring
            allConv = allConv -  BLmean;
            allConv = allConv ./ BLstd;
            
            allConvTun  = allConv( tunTrials,:); % extract only tuned trials
            %             allConvNun  = allConv(~tunImg,:);
            
            %             allTunSU    = [allTunSU, {allConvTun}];
            %             allNunSU    = [allNunSU, {allConvNun}];
            allTun      = [allTun; allConvTun];
            
            
        end
        save('quirogaTuning.mat', 'image_cell');
    end
end

save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\declan tuning\allTun_quiroga_blPreCue.mat', 'allTun');
end % END OF FUNCTION

clf
hold on;
for tt = 1:32
    plot(allTun(tt,:), 'linew', 2, 'color', [0 0 0 0.5]);
end

