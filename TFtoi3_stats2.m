%% Here I z-score first and then take the difference
%% WITHIN TRIAL
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\TF_preCue.mat', 'allSpks');
clearvars -except allSpks

% PREALLOCATE
% LOW
diffTFlow_all = [];
idxNumLow_all = [];
powLow_all    = [];

% HIGH
diffTFhigh_all = [];
idxNumHigh_all = [];
powHigh_all    = [];

subsRel = 0; % 0 for subs, 1 for rel

for su = 1 : length(allSpks)
    lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', su, length(allSpks), su/length(allSpks)*100);
    
    idxTrlLw = allSpks(su).idxTrlSingLw;
    idxTrlHi = allSpks(su).idxTrlSingHi;
    
    %% LOW
    if ~any(isnan(allSpks(su).idxTrlSingLw))
        if sum(allSpks(su).idxTrlSingLw) > 0
            
            % SAVE STUFF FOR PERMUTATION TEST
            idxNum    = sum(allSpks(su).idxTrlSingLw);
            powLow    = cat(3, allSpks(su).allPowLow(:,:,idxTrlLw), allSpks(su).allPowLow(:,:, ~idxTrlLw));
            
            %% zscore with all trials as BL
            powLowBase = squeeze(nanmean(powLow,2)); % first collapse over time
            baseMean   = mean(powLowBase,2);  %mean over trials; this is the mean for each frequency
            baseStd    = std(powLowBase,0,2);
            
            % ZSCORING
            for trl = 1:size(powLow,3)
                powLow(:,:,trl) = (powLow(:,:,trl)-baseMean)./baseStd;
            end
            
            powLowIdx = powLow(:,:,1:idxNum);
            powLowNdx = powLow(:,:,idxNum+1:end);
            
            %% TAKE DIFFERENCE
            % TF of difference (ABSOLUTE OR RELATIVE)
            switch subsRel
                case 0
                    diffTF = squeeze(nanmean(powLowIdx,3)) - squeeze(nanmean(powLowNdx,3));
                case 1
                    diffTF = squeeze(nanmean(powLowIdx,3)) ./ squeeze(nanmean(powLowNdx,3));
            end
                        
            
            diffTFlow_all = [diffTFlow_all {diffTF}];
            idxNumLow_all = [idxNumLow_all, idxNum];
            powLow_all    = [powLow_all, {powLow}];
            
        end
    end
    
   %% HIGH
    if ~any(isnan(allSpks(su).idxTrlSingHi))
        if sum(allSpks(su).idxTrlSingHi) > 0
            
            % SAVE STUFF FOR PERMUTATION TEST
            idxNum    = sum(allSpks(su).idxTrlSingHi);
            powHigh    = cat(3, allSpks(su).allPowHigh(:,:,idxTrlHi), allSpks(su).allPowHigh(:,:, ~idxTrlHi));

                       %% zscore with all trials as BL
            powHighBase = squeeze(nanmean(powHigh,2)); % first collapse over time
            baseMean   = mean(powHighBase,2);  %mean over trials; this is the mean for each frequency
            baseStd    = std(powHighBase,0,2);
            
            % ZSCORING
            for trl = 1:size(powHigh,3)
                powHigh(:,:,trl) = (powHigh(:,:,trl)-baseMean)./baseStd;
            end
            
            powHighIdx = powHigh(:,:,1:idxNum);
            powHighNdx = powHigh(:,:,idxNum+1:end);
            
            %% TAKE DIFFERENCE
            % TF of difference (ABSOLUTE OR RELATIVE)
            switch subsRel
                case 0
                    diffTF = squeeze(nanmean(powHighIdx,3)) - squeeze(nanmean(powHighNdx,3));
                case 1
                    diffTF = squeeze(nanmean(powHighIdx,3)) ./ squeeze(nanmean(powHighNdx,3));
            end
            
            diffTFhigh_all = [diffTFhigh_all, {diffTF}];
            idxNumHigh_all = [idxNumHigh_all, idxNum];
            powHigh_all    = [powHigh_all, {powHigh}];        
        end
    end
    
    
    fprintf(repmat('\b',1,lineLength))
end

%% EMPIRICAL DIFFERENCE 
% LOW
diffTFlow_all = cat(3,diffTFlow_all{:});
diffTFlow_all = squeeze(nanmean(diffTFlow_all,3));

% HIGH
diffTFhigh_all = cat(3,diffTFhigh_all{:});
diffTFhigh_all = squeeze(nanmean(diffTFhigh_all,3));

%% PERMUTATION TEST: LOW
nperm = 100;
powLowPerm = zeros(nperm, 45, 1000);
for perm = 1:nperm
        lineLength = fprintf('%d permutations out of %d permutations done (%.2f%%).\n', perm, nperm, perm/nperm*100);

    diffTFperm_comp = [];
    
    for comp = 1 : size(powLow_all,2)
        
        numIdx = idxNumLow_all(comp);
        seshTF = powLow_all{comp}(:,:,randperm(size(powLow_all{comp},3))); %% SHUFFLE
             
        %% USE ABSOLUTE OR RELATIVE DIFFERENCE
        switch subsRel
            case 0
                diffTFperm = squeeze(nanmean(seshTF(:,:,1:numIdx),3) - nanmean(seshTF(:,:, numIdx+1:end),3));
            case 1
                diffTFperm = squeeze(nanmean(seshTF(:,:,1:numIdx),3)) ./ squeeze(mammean(seshTF(:,:, numIdx+1:end),3));
        end
        
        diffTFperm_comp(:,:,comp) = diffTFperm;        
    end
    
    diffTFperm_comp      = squeeze(nanmean(diffTFperm_comp,3));
    powLowPerm(perm,:,:) = diffTFperm_comp;
end


% mean_h0 = squeeze(nanmean(powLowPerm,1));
% std_h0  = squeeze(nanstd(powLowPerm,1));
zval    = abs(norminv(0.05));

max_vals  = zeros(nperm,2);
bigIsland = zeros(nperm,1);
for perm = 1:nperm

    %% THE PERMUTED TFs ARE ALREADY Z-NORMALIZED. WHY DO WE NORMALIZE HERE AGAIN??
%     permMap = squeeze(powLowPerm(perm,:,:));
%     permMap = (permMap-mean_h0)./std_h0;
    
    permMap = squeeze(powLowPerm(perm,:,:)); % WHY NOT THIS INSTEAD?

    %% FIND MAX VALUE
    temp = sort( reshape(permMap,1,[] ));
    temp(isnan(temp)) = [];
    max_vals(perm,:) = temp([1 end]);
    
    %% FIND ISLANDS
    permMap(abs(permMap)<zval) = 0;
    islands = bwconncomp(permMap);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        bigIsland(perm) = max(tempclustsizes);
    end
end

% find the lower threshold bound based on percentage
tmp = sort(max_vals(:,1));
thresh_lo = tmp( round(0.05*nperm) ); % SAME AS: prctile(tmp, 5)

% repeat for upper threshold bound
tmp = sort(max_vals(:,2));
thresh_hi = tmp( round((1-0.05)*nperm) );

%%
figure(1); clf; hold on;
contmap = powLowDiff>=thresh_hi | powLowDiff <= thresh_lo;
imagesc(powLowDiff)

% threshold real data
diffmap = squeeze(nanmean(powLowEmp,3));
zmap = squeeze(nanmean(powLowEmp,3));
zmap(zmap>thresh_lo & zmap<thresh_hi) = 0;


%% PERMUTATION TEST: HIGH
nperm = 100;
powHighPerm = zeros(nperm, 81, 1000);
for perm = 1:nperm
        lineLength = fprintf('%d permutations out of %d permutations done (%.2f%%).\n', perm, nperm, perm/nperm*100);

    diffTFperm_comp = [];
    
    for comp = 1 : size(powHigh_all,2)
        
        numIdx = idxNumHigh_all(comp);
        seshTF = powHigh_all{comp}(:,:,randperm(size(powHigh_all{comp},3))); %% SHUFFLE
             
        %% USE ABSOLUTE OR RELATIVE DIFFERENCE
        switch subsRel
            case 0
                diffTFperm = squeeze(nanmean(seshTF(:,:,1:numIdx),3) - nanmean(seshTF(:,:, numIdx+1:end),3));
            case 1
                diffTFperm = squeeze(nanmean(seshTF(:,:,1:numIdx),3)) ./ squeeze(mammean(seshTF(:,:, numIdx+1:end),3));
        end
        
        % z transform
        for freq = 1:size(diffTFperm,1)
            avg                = nanmean(diffTFperm(freq,:),2);
            mStdev             = nanstd(diffTFperm(freq,:),0,2);
            diffTFperm(freq,:) = (diffTFperm(freq,:)-avg)/mStdev;
        end
        
        diffTFperm_comp(:,:,comp) = diffTFperm;        
    end
    
    diffTFperm_comp      = squeeze(nanmean(diffTFperm_comp,3));
    powHighPerm(perm,:,:) = diffTFperm_comp;
end

zval    = abs(norminv(0.05));

max_vals  = zeros(nperm,2);
bigIsland = zeros(nperm,1);
for perm = 1:nperm
    
    permMap = squeeze(powHighPerm(perm,:,:)); 

    %% FIND MAX VALUE
    temp = sort( reshape(permMap,1,[] ));
    temp(isnan(temp)) = [];
    max_vals(perm,:) = temp([1 end]);
    
    %% FIND ISLANDS
    permMap(abs(permMap)<zval) = 0;
    islands = bwconncomp(permMap);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        bigIsland(perm) = max(tempclustsizes);
    end
end

% LOWER THRESHOLD
tmp = sort(max_vals(:,1));
thresh_lo = tmp( round(0.05*nperm) ); % SAME AS: prctile(tmp, 5)

% repeat for upper threshold bound
tmp = sort(max_vals(:,2));
thresh_hi = tmp( round((1-0.05)*nperm) );