clear

folderpath = '\\analyse4.psy.gla.ac.uk\project0309\Luca\official code'; % wherever you have downloaded this folder
load([folderpath, '\data\allSpks.mat'], 'allSpks');
pPermNorm = [];
pPerm     = [];
nperm     = 1000;

% GAUSSIAN KERNEL
mlength = 251;
mKernel = gausswin(mlength);

tic
for su = 1:length(allSpks)
    timePassed = toc;
    ETA = timePassed/(su-1)*(length(allSpks)-su);
    lineLength = fprintf('Analysing SU#%d | Time passed: %.1f min | Done in: %.1f min \n', su, (timePassed/60), (ETA/60));
    
    
    spks    = round(allSpks(su).spks);
    
    
    %% HITS AND MISSES + WHOLE TRIAL
    encTrig = round(allSpks(su).encTrigger(:,[1 3])*1000);
    retTrig = round(allSpks(su).retTrigger(:,[1 3])*1000);
    
    enc = encTrig(:,1);
    ret = retTrig(:,1);
    
    block       = 1;
    encBlock    = [];
    for trl = 1:size(enc,1)-1
        encTrl = enc(trl);
        encNxt = enc(trl+1);
        
        retTrls = abs(encNxt - ret);
        
        encBlock(trl,1) = block;
        
        if any(encNxt-encTrl > retTrls)
            block = block + 1;
        end
        
    end
    encBlock(end+1) = block;
    encTrig = [encTrig, encBlock];
    
    nBlocks = max(encTrig(:,3));
    
    if nBlocks == 1 % there is only one block
        allSpks(su).pTCnorm = NaN;
        allSpks(su).pTC     = NaN;
        continue
    end
    
    for block = 1:nBlocks
        blockStartEnd(block,:) = [min(encTrig( encTrig(:,3)==block ,1)), max(encTrig( encTrig(:,3)==block ,2))]; % start and end of the block
    end
    
    spksBlock = {};
    for block = 1:nBlocks
        spksBlock{block} = spks( spks >= blockStartEnd(block,1) & spks <= blockStartEnd(block,2) ) - blockStartEnd(block,1) + 1;
    end
    
    % DETERMINE LONGEST BLOCK LENGTH
    blength = blockStartEnd(:,2)-blockStartEnd(:,1);
    blengthMax = max(blength);
    
    %% No Normalisation
    %  SEGMENT BLOCKS INTO 40 EQUALLY SIZED BINS
    spkConv = [];
    for block = 1 : nBlocks
        
        % bin spikes for each block and perform gaussian convolution
        dt2                  = [0:5:blockStartEnd(block,2)- blockStartEnd(block,1)];
        [x2,~]               = histcounts(spksBlock{block}, dt2);
        spkConv{block}       = conv(x2,mKernel,'same');
        
        dtNorm               = round(linspace(1,size(spkConv{block},2),41));
        if dtNorm(end) > size(spkConv{block},2)
            dtNorm(end) = size(spkConv{block},2);
        end
        
        dat = spkConv{block};
        for mm = 1:length(dtNorm)-1
            snipDat = mean(dat(dtNorm(mm):dtNorm(mm+1)));
            spkConvNorm(block, mm) = snipDat;
        end
        
    end
    
    % KRUSKALWALLIS
    [pTC,~,~] = kruskalwallis(spkConvNorm, [], 'off');
    allSpks(su).pTCnorm = pTC;
    
    % PERMUTATIONTEST (HOW MANY TC UNDER H0 -> pPermNorm)
    parfor perm = 1 : nperm
        spkConvNormPerm = [];
        for block = 1 : nBlocks
            
            dtNorm = round(linspace(1,size(spkConv{block},2),41));
            if dtNorm(end) > size(spkConv{block},2)
                dtNorm(end) = size(spkConv{block},2);
            end
            
            dat2Shift = spkConv{block};
            dat2Shift = circshift(dat2Shift,randi([1,length(dat2Shift)],1),2);
            for mm = 1:length(dtNorm)-1
                snipDat = mean(dat2Shift(dtNorm(mm):dtNorm(mm+1)));
                spkConvNormPerm(block, mm) = snipDat;
            end
            
        end
        tt = kruskalwallis(spkConvNormPerm, [], 'off');
        pPermNorm(su,perm) = tt;
        
    end
    
    %% Normalization
    % SEGMENT BLOCKS INTO 40 BINS WITHOUT NORMALIZING TO THE SAME LENGTH
    spkConvNorm = [];  % <. that shouldnt be nonorm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% and also adjust for skipped
    spkConv     = [];
    dtNorm = round(linspace(1,blengthMax,41));
    for block = 1 : nBlocks
        % bin spikes for each block and perform gaussian convolution
        dt2              = [0:1:blockStartEnd(block,2)- blockStartEnd(block,1)];
        [x2,~]           = histcounts(spksBlock{block}, dt2);
        dat              = conv(x2,mKernel,'same');
        dat              = [dat, NaN( 1, blengthMax - size(dat,2) )];
        spkConv(block,:) = dat;
        
        for mm = 1:length(dtNorm)-1
            snipDat                = mean(dat(dtNorm(mm):dtNorm(mm+1)));
            spkConvNorm(block, mm) = snipDat;
        end
    end
    
    % KRUSKALWALLIS
    [pTC,~,~]           = kruskalwallis(spkConvNorm, [], 'off');
    allSpks(su).pTC     = pTC;
%     allSpks(su).spkConv = spkConvNorm;
    
    % PERMUTATIONTEST (HOW MANY TC UNDER H0 -> pPermNorm)
    parfor perm = 1 : nperm
        spkConvNormPerm = [];
        for block = 1 : nBlocks
            
            % circular shift current block
            dat2Shift                   = spkConv(block,:);
            numNan                      = sum(isnan(dat2Shift)); % how many nans does this block have?
            dat2Shift(isnan(dat2Shift)) = []; % get rid of nans before circular shuffling
            dat2Shift                   = circshift(dat2Shift,randi([1,length(dat2Shift)],1),2); % circular shuffle
            dat2Shift                   = [dat2Shift, NaN(1,numNan)]; % add nans again
            
            % bin into 40 bins
            for mm = 1:length(dtNorm)-1
                snipDat                    = mean(dat2Shift(dtNorm(mm):dtNorm(mm+1)));
                spkConvNormPerm(block, mm) = snipDat;
            end
            
        end
        pPerm(su,perm) = kruskalwallis(spkConvNormPerm, [], 'off');
        
    end
    fprintf(repmat('\b', 1, lineLength));
end

save([folderpath, '\TC-KW_norm_noNorm.mat'], 'pPermNorm', 'pPerm', 'allSpks')



%% Level 2 permutation test
for norm=1:2
    
    switch norm
        case 1
            permDist  = pPerm;     % NO NORMALISATION
            pTC = [allSpks.pTC];
        case 2
            permDist  = pPermNorm; % WITH NORMALISATION
            pTC = [allSpks.pTCnorm];
    end
    
    permDist = permDist<0.05;
    
    numPerm = zeros(1,nperm);
    for perm = 1:nperm % repeat nperm-times
        for su = 1:size(permDist,1) % randomly choose a permutation value for each su
            randIdx = randperm(size(permDist,2),1);
            numPerm(perm) = numPerm(perm) + permDist(su, randIdx);
        end
    end
    
    %% VISUALISATION
    idx    = pTC <= 0.05; % which SU are classified as time cells?
    allESN = [allSpks.ESN];
    allESN(allESN == 2 | allESN == 1) = 1;
    sum(allESN(idx)) % how many TC are ESN?
    
    figure(1); clf; hold on;
    histogram(numPerm, 'normalization', 'probability')
    xlabel('Number of TC under the H_0')
    ylabel('Probability')
    plvl = mean(numPerm>=sum(pTC <=0.05));
    title(plvl);
    plot([sum(pTC <= 0.05) sum(pTC <= 0.05)], get(gca, 'YLim'), '--', 'color', [1 0 0], 'linew', 2)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',35,'FontWeight','Bold', 'LineWidth', 2);
    
    
    figure(2); clf;
    histogram(pTC, [0:0.01:1], 'normalization', 'probability')
    hold on;
    plot([0.05 0.05], get(gca, 'YLim'), '--', 'color', [1 0 0], 'linew', 2)
    xlabel('p values')
    ylabel('proportion')
    title(sprintf('How many SU (#625) are TC (#%d | %d are pIU)?', sum(pTC <= 0.05), sum(allESN(pTC <= 0.05))))
    mAx = gca;
    mAx.FontWeight = 'bold';
    mAx.FontSize = 20;
    
    %% ESN AND TC INTERACTION
    numTC = sum(pTC <= 0.05); % number of TC XX put above
    
    % How many ESNs are expected to be drawn from numTC number of drawings under the H0?
    for perm = 1:nperm
        
        allESNperm    = allESN(randperm(size(allESN,2)));
        permNum(perm) = sum(allESNperm(1:numTC));
        
    end
    
    numESN = sum(allESN(pTC <= 0.05)); % number of TC that are also ESN XX put above
    pVal = mean(permNum >= numESN)
    
    
end
