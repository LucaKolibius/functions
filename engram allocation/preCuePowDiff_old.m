clearvars -except trigALL
lfpDir         = dir('X:\George\Analysis\Data Continuous\sub-*.mat');
idxPow2        = [];
ndxPow2        = [];
idxPow         = [];
ndxPow         = [];
dxdiff         = [];
trlChanPowIdx  = [];
trlChanPowNdx  = [];
ndxPowGED      = [];  % weighted component
idxPowGED      = [];
ndxPowGEDperm  = [];
idxPowGEDperm  = [];
nperm          = 1000;


for ses = 1 : size(lfpDir,1)
    disp(ses)
    
    % get current bidsID & sesh
    bidsID = lfpDir(ses).name; sesh = bidsID;
    bidsID = bidsID(1:8);
    sesh   = sesh(10:12);
    
    if ~strcmp(sesh(3), 'b')
        sesh(3) = [];
    end
    
    % load the data
    load([lfpDir(ses).folder, filesep, lfpDir(ses).name], 'data_micro');
    
    % select the LFP (ignore spikes / trigger / hits)
    labs   = data_micro.label;
    noSpks = cellfun(@(x) regexp(x, 'SPKS'), labs, 'un', 0);
    noTrig = cellfun(@(x) regexp(x, 'TRIG'), labs, 'un', 0);
    noHits = cellfun(@(x) regexp(x, 'enc+'), labs, 'un', 0);
    lfpNam =  cellfun(@isempty, noSpks);
    noTrig = ~cellfun(@isempty, noTrig);
    noHits = ~cellfun(@isempty, noHits);
    
    lfpNam = lfpNam - noTrig - noHits;                          % index all the names of LFPs
    LFPlab = labs(logical(lfpNam));                             % names of all the LFPs (MB1, MB2, MB3...)
    bndLab = unique(cellfun(@(x) x(1:end-1), LFPlab, 'un', 0)); % only the unique bundles
    
    for bnd = 1 : size(bndLab,1)
        
        % which rows of lfp data do we need from data_micro?
        microIdx = cellfun(@(x) regexp(x,bndLab(bnd)), {labs}, 'un', 0);
        microIdx = microIdx{1};
        microIdx = cellfun(@(x) find(x == 1), microIdx, 'un', 0);
        microIdx = ~cellfun(@isempty, microIdx);
        
        % get the LFP from data_micro
        dat = data_micro.trial{1}(microIdx,:);
        
        % FIND THE TRIALS OF THAT BUNDLE THAT ARE INDEXING (IF ANY)
        % extract all SU that fire during that session + bundle
        allBund  = {trigALL.wirename};
        allBund  = cellfun(@(x) x(1:end-1), allBund, 'un', 0);
        sameBund = and(and(contains({trigALL.bidsID}, bidsID), strcmp({trigALL.sesh}, sesh) ), contains(allBund, bndLab(bnd))); % same subject + session + bundle
        
        idxTrl = any(vertcat(trigALL(sameBund).idxTrl),1); % these are the trials that are indexed in that wire
        
        % SOMETIMES THERE ARE NO SPIKES. THEN I JUST TAKE THE BEST MW FROM
        % THAT SESSION AND PARTICIPANT
        sameBidsSesh = and( contains({trigALL.bidsID}, bidsID), contains({trigALL.sesh}, sesh));
        sameBidsSesh = find(sameBidsSesh == 1);
        if isempty(sameBidsSesh)
            continue
        end
        
        sameBidsSesh = sameBidsSesh(1);
        
        % if the session does not have any spikes, no trials are indexed
        % trials. But I need to get the number of trials from a different
        % bundle in that same session
        if isempty(idxTrl)
            idxTrl = logical(zeros(1,size(trigALL(sameBidsSesh).idxTrl,2)));
        end
        
        encTrig    = round(trigALL(sameBidsSesh).encTrigger(trigALL(sameBidsSesh).hitsIdx) * 1000); % encoding trigger
        encTrig    = [encTrig-1000 encTrig];
        
        hz         = linspace(0,1000, 1001);
        trlPow     = zeros(size(idxTrl,2), 1001);
        covAll     = [];
        trlChanPow = [];
        % IMPORANT: I USE THE AVERAGE PWSPCTRM OVER ALL MW. THIS IS NOT
        % OPTIMAL. VISUALIZE.
        for trl = 1 : size(idxTrl, 2)
            datSnip             = dat(:, encTrig(trl,1) : encTrig(trl,2));
            trlPow(trl,:)       = mean(abs(fft(datSnip)),1);
            trlChanPow(trl,:,:) = abs(fft(datSnip));
            
            datSnip         = bsxfun(@minus,datSnip,mean(datSnip,2)); % demean
            covAll(trl,:,:) = datSnip*datSnip';
        end
        
        %% GED IDEA
        if sum(idxTrl) > 0
            
            numIdx = sum(idxTrl);                                           % number of indexed trials
            
            for perm = 1 : nperm
                shuf   = randperm(size(covAll,1));                          % shuffling index
                tmpCov = covAll(shuf, :,:);                                 % shuffle all covariance matrices
                
                covPermIdx = squeeze(mean(tmpCov(1:numIdx,     :, :),1));   % generate averaged, permuted covariance matrix for indexed trials
                covPermNdx = squeeze(mean(tmpCov(numIdx+1:end, :, :),1));   % generate averaged, permuted covariance matrix for nondxed trials
                
                [evecP, evalP] = eig(covPermIdx, covPermNdx);               % GED
                [~,maxcomp]    = sort(diag(evalP));                         % index of the biggest component
                evec           = evecP(:,maxcomp(end));
                
                for trl = 1 : numIdx
                    
                    datSnip      = dat(:, encTrig(shuf(trl),1) : encTrig(shuf(trl),2));
                    datSnip      = (datSnip' * evec)';                          % weight the MW according to eigenvalues
                    
                    temp(trl,:)  = abs(fft(datSnip));
                    
                end
                idxPowGEDperm    = [idxPowGEDperm; mean(temp,1)];
                
            end
            
            
            covIdx = squeeze(mean(covAll( idxTrl,:,:),1));
            covNdx = squeeze(mean(covAll(~idxTrl,:,:),1));
            
            [evecs,evals] = eig(covIdx,covNdx); % compute eigenvalues
            
            % find best component and compute filter projection
            [~,maxcomp] = sort(diag(evals)); % index of the biggest component
            evec        = evecs(:,maxcomp(end));
            
            for trl = 1 : size(idxTrl, 2)
                datSnip  = dat(:, encTrig(trl,1) : encTrig(trl,2));
                
                %                 % only take the MW with the highest eigenvalue
                %                 [~,pos] = max(evec);
                %                 datSnip2 = datSnip(pos,:);
                
                % weight the MW according to eigenvalues
                datSnip  = (datSnip' * evec)';
                
                switch idxTrl(trl)
                    case 0
                        ndxPowGED  = [ndxPowGED;  abs(fft(datSnip))];   % weighted component
                        %                         ndxPowGED2 = [ndxPowGED2; abs(fft(datSnip2))];  % MW with highest eval
                    case 1
                        idxPowGED  = [idxPowGED;  abs(fft(datSnip))];
                        %                         idxPowGED2 = [idxPowGED2; abs(fft(datSnip2))];
                end
            end
        end
        %%
        % use every indexed trial and non-indexed trial as a random effect
        idxPow = [idxPow; trlPow( idxTrl,:)];
        ndxPow = [ndxPow; trlPow(~idxTrl,:)];
        
        % log the average difference per bundle in indexed vs. non-indexed
        % trials
        dxdiff  = [dxdiff;   nanmean(trlPow( idxTrl,:) ,1) - nanmean(trlPow(~idxTrl,:) ,1)];
        
        % log the mean power response for indexed and non indexed trials in
        % that bundle
        idxPow2 = [idxPow2, {nanmean(trlPow( idxTrl,:) ,1)} ];
        ndxPow2 = [ndxPow2, {nanmean(trlPow(~idxTrl,:) ,1)} ];
        
        %%
        trlChanPowIdx = [trlChanPowIdx; trlChanPow( idxTrl,:,:)];
        trlChanPowNdx = [trlChanPowNdx; trlChanPow(~idxTrl,:,:)];
    end
end


%%
for trl = 1:5:50 % size(loggItIDX,1)
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for chan = 1:8
        subplot(8,1,chan)
        plot(squeeze(trlChanPowIdx(trl, chan, :)), 'linew', 2);
        title('IDX')
        xlim([0 200])
    end
    
    
    drawnow
end

%%
figure('units','normalized','outerposition',[0 0 1 1])
hold on
for ii = 1:200
    plot(hz, idxPow(ii,:))
end

newIdx = idxPow;
newIdx = newIdx./max(newIdx,[], 2);

newNdx = ndxPow;
newNdx = newNdx./max(newNdx,[], 2);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
plot(hz, nanmean(newIdx,1), 'color', 'b', 'linew', 3);
plot(hz, nanmean(newIdx,1) - 1.96 * nanstd(newIdx,1), 'color', 'b', 'linew', 3);
plot(hz, nanmean(newIdx,1) + 1.96 * nanstd(newIdx,1), 'color', 'b', 'linew', 3);

plot(hz, nanmean(newNdx,1), 'color', 'r', 'linew', 3);
plot(hz, nanmean(newNdx,1) - 1.96 * nanstd(newNdx,1), 'color', 'r', 'linew', 3);
plot(hz, nanmean(newNdx,1) + 1.96 * nanstd(newNdx,1), 'color', 'r', 'linew', 3);

xlim([0 100])

%% GED VISU
% figure('units','normalized','outerposition',[0 0 1 1]); hold on;
% hz = linspace(0,1000, 1001);
% idxPowGEDnorm = idxPowGED ./ max(idxPowGED, [], 2);
% subplot(211)
% imagesc(idxPowGEDnorm)
% xlim([0 200])
%
% subplot(212)
% ndxPowGEDnorm = ndxPowGED ./ max(ndxPowGED, [], 2);
% imagesc(ndxPowGEDnorm)
% xlim([0 200])

% figure('units','normalized','outerposition',[0 0 1 1]); hold on;
% idxPowGED2norm = idxPowGED2% ./ max(idxPowGED2, [], 2);
% subplot(211)
% ii = imagesc(idxPowGED2norm)
% yy = colorbar;
% caxis manual
% % caxis([0 0.5])
% caxis([0 3000])
% xlim([0 200])
%
% subplot(212)
% ndxPowGED2norm = ndxPowGED2% ./ max(ndxPowGED2, [], 2);
% imagesc(ndxPowGED2norm)
% caxis manual
% % caxis([0 0.5])
% caxis([0 3000])

% xlim([0 200])

% weighted component
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
plot(hz, mean(idxPowGED,1), 'linew', 3);
plot(hz, mean(ndxPowGED,1), 'linew', 3);

thUP = prctile(idxPowGEDperm, 95, 1);
plot(hz, thUP, 'linew', 2, 'color', 'r');
xlim([0 100])

%% artefacts
% should be implemented using lfp2microLFP

%% permtest
nperm = 10000;
allDat = [newIdx; newNdx];
permDx = [];
for perm = 1:nperm
    randdx = randperm(size(allDat,1));
    
    permDx(perm,:) = mean(allDat(randdx(1:size(idxPow,1)),:),1);
end

upTH = prctile(permDx,95,1);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
plot(hz, nanmean(newIdx,1), 'color', 'b', 'linew', 3);
plot(hz, upTH, 'color', 'r', 'linew', 3);