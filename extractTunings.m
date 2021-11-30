%% extracts spiketimes from CN in declans data
function extractTunings
cd('\\analyse4.psy.gla.ac.uk\project0309\Declans_Tuning_Data\tuning');
allSbj = dir('sub-00*');

allTun = [];
skipped  = [];

allTunSU    = [];
allNunSU    = [];

%% gaussian kernel
% mlength = [-0.075:0.0002:0.075];
% mSigma  = 0.02; % 20ms
% mKernel = normpdf(mlength,0,mSigma);
% mKernel = mKernel/max(mKernel); % normalize peak to 1
mlength = 251;
mKernel = gausswin(mlength);

dt      = -0.625:0.001:1.125+0.001;
dt      = dt-0.0005; % center around 0


for sbj = 1:size(allSbj,1)
    cd([allSbj(sbj).folder, filesep, allSbj(sbj).name]);
    
    allSesh = dir('20*');
    
    for sesh = 1:size(allSesh,1)
       allSU = {};
       
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
        numSU    = size(resp_mat,1);
        
        for su = 1:numSU
            allSU(su,:) = insertSpiketimes(ttl_vec, allCl{su,2}, 1, [-0.625 1.125]);
        end
        
        for su = 1:numSU
            idxStim  = find(resp_mat(su,:) == 1); % resp_mat is SUxStimulus. find which stimuli for which SU are tuned
            
            if isempty(idxStim)
                continue
            end
            
            idxTrl                          = zeros(max([image_cell{:,2}]),1);
            idxTrl([image_cell{idxStim,2}]) = 1; % based on the tuned stimuli, which trials are tuned?
            idxTrl                          = logical(idxTrl);
            
            allConvRaw = [];
            for trl = 1:size(allSU,2)
                
                [x,~]   = histcounts(allSU{su,trl}, dt);
                spkConv = conv(x, mKernel, 'same');
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
            
            allConvTun  = allConv( idxTrl,:); % extract only tuned trials
            allConvNun  = allConv(~idxTrl,:);
            
%             figure(2); clf; %hold on;
%             for tun = 1:6
%                 subplot(2,3,tun)
%                 plot(allConvTun(tun,:), 'color', [0 0 0 0.5], 'linew', 2); ylim([-2 7])
% %                 pause(1.5);
%             end
%             sgtitle('Normalized Tuned');
            
%             %% raw
%             allConvRaw = allConvRaw(idxTrl,627:end-125);
%              figure(1); clf; %hold on;
%             for tun = 1:6
%                 subplot(2,3,tun)
%                 plot(allConvRaw(tun,:), 'color', [0 0 0 0.5], 'linew', 2); ylim([-2 7])
% %                 pause(1.5);
%             end
            
%             maxVals = max(allConvTun,[],2);
%             figure(3);
%             histogram(maxVals, 0:0.5:8)
            
            allTun      = [allTun; allConvTun]; % don't mean(allConvTun,1)
            
            allTunSU    = [allTunSU, {allConvTun}];
            allNunSU    = [allNunSU, {allConvNun}];
            
        end % SU END
        
    end % SESSION END
end % SUBJECT END

save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\declan tuning\allTun_blPreCue.mat', 'allTun');
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\declan tuning\allTunNun_SU.mat' , 'allTunSU', 'allNunSU');
end % END OF FUNCTION
