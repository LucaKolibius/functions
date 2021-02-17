%% extracts spiketimes from CN in declans data
function extractTunings
cd('\\analyse4.psy.gla.ac.uk\project0309\Declans_Tuning_Data\tuning');
allSbj = dir('sub-00*');

allTun = [];
skipped  = [];

%% gaussian kernel
mlength = [-0.075:0.0002:0.075];
mSigma  = 0.02; % 20ms
mKernel = normpdf(mlength,0,mSigma);
mKernel = mKernel/max(mKernel); % normalize peak to 1
dt      = -0.375:0.001:1.375+0.001;
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
        load('allrespData', 'resp_mat');
        catch
            skipped = [skipped; sbj, sesh];
            continue
        end
        numSU    = size(resp_mat,1);
        
        for su = 1:numSU
            allSU(su,:) = insertSpiketimes(ttl_vec, allCl{su,2}, 1, [-0.375 1.375]);
        end
        
        for su = 1:numSU
            idxStim  = find(resp_mat(su,:) == 1); % resp_mat is SUxStimulus. find which stimuli for which SU are tuned
            
            if isempty(idxStim)
                continue
            end
            
            idxTrl   = [image_cell{idxStim,2}]; % based on the tuned stimuli, which trials are tuned?
            
            allConv = [];
            for trl = 1:size(allSU,2)
                
                [x,~]                = histcounts(allSU{su,trl}, dt);
                spkConv              = conv(mKernel, x);
                
                spkConv(1:750)       = [];
                spkConv(end-749:end) = [];
                
                allConv = [allConv; spkConv];
                
            end
            BL = mean(allConv,2); % over time
            BLmean = mean(BL);
            BLstd  = std(BL, 0);
            
            allConv = allConv -  BLmean;
            allConv = allConv ./ BLstd;
            
            allConvTun  = allConv( idxTrl,:);
            
            allTun      = [allTun; allConvTun]; % don't mean(allConvTun,1)
            
        end % SU END
        
    end % SESSION END
end % SUBJECT END

save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\declan tuning\allTun.mat', 'allTun');
end % END OF FUNCTION
