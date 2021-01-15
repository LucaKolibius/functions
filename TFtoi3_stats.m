%% WITHIN TRIAL
% load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\TFtoi3.mat', 'allSpks');
clearvars -except allSpks
powLowDiff    = [];
powLowComp    = [];
powLowCompIdx = [];

powHighDiff    = [];
powHiComp    = [];
powHiCompIdx = [];

for su = 1 : length(allSpks)
    lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', su, length(allSpks), su/length(allSpks)*100);
    
    idxTrlLw = allSpks(su).idxTrlSingLw;
    idxTrlHi = allSpks(su).idxTrlSingHi;
    
    %% LOW
    if ~any(isnan(allSpks(su).idxTrlSingLw))
        if sum(allSpks(su).idxTrlSingLw) > 0
            
            %% SAVE STUFF FOR PERMUTATION TEST
            idxNum    = sum(allSpks(su).idxTrlSingLw);
            powLowAll = cat(3, allSpks(su).allPowLow(:,:,idxTrlLw), allSpks(su).allPowLow(:,:, ~idxTrlLw));


            
%             %% NORMALIZE IDX-TF AND NDX-TF AND THEN TAKE DIFFERENCE
%             powLowIdx = squeeze(nanmean(allSpks(su).allPowLow(:,:,idxTrlLw),3));
%             powLowNdx = squeeze(nanmean( allSpks(su).allPowLow(:,:, ~idxTrlLw),3));
%             for freq = 1:size(powLowIdx,1)
%                 powLowIdx(freq,:) = 10*log10(powLowIdx(freq,:)./squeeze(nanmean(powLowIdx(freq,:),2)));
%             end
%             
%             for freq = 1:size(powLowNdx,1)
%                 powLowNdx(freq,:) = 10*log10(powLowNdx(freq,:)./squeeze(nanmean(powLowNdx(freq,:),2)));
%             end
%             
%             powLowDiff = [powLowDiff, {powLowIdx-powLowNdx}];
            
            %%
            
            %% NORMALIZE EVERY TRIAL (PRODUCES NO RESULT)
            
            % DB NORMALIZE ALL TRIALS
            for trl = 1:size(powLowAll,3)
                for freq = 1:size(powLowAll,1)
                    powLowAll(freq,:, trl) = 10*log10(powLowAll(freq,:, trl)./squeeze(nanmean(powLowAll(freq,:, trl),2)));
                end
            end
            
            powLowComp    = [powLowComp {powLowAll}];
            powLowCompIdx = [powLowCompIdx, idxNum];

            powLowIdx  = powLowAll(:,:,1:idxNum);
            powLowNdx  = powLowAll(:,:,idxNum+1:end);
            powLowDiff = [powLowDiff, {squeeze(nanmean(powLowIdx,3))-squeeze(nanmean(powLowNdx,3))} ]; % should be faster than concatenating to third dimension
                        
        end
    end
    
   %% HIGH
    if ~any(isnan(allSpks(su).idxTrlSingHi))
        if sum(allSpks(su).idxTrlSingHi) > 0
            
            %% SAVE STUFF FOR PERMUTATION TEST
            idxNum    = sum(allSpks(su).idxTrlSingHi);
            powHighAll = cat(3, allSpks(su).allPowHigh(:,:,idxTrlHi), allSpks(su).allPowHigh(:,:, ~idxTrlHi));

            %% NORMALIZE EVERY TRIAL
            % DB NORMALIZE ALL TRIALS
            for trl = 1:size(powHighAll,3)
                for freq = 1:size(powHighAll,1)
                    powHighAll(freq,:, trl) = 10*log10(powHighAll(freq,:, trl)./squeeze(nanmean(powHighAll(freq,:, trl),2)));
                end
            end
            
            powHiComp    = [powHiComp {powHighAll}];
            powHiCompIdx = [powHiCompIdx, idxNum];

            powHighIdx  = powHighAll(:,:,1:idxNum);
            powHighNdx  = powHighAll(:,:,idxNum+1:end);
            powHighDiff = [powHighDiff, {squeeze(nanmean(powHighIdx,3))-squeeze(nanmean(powHighNdx,3))} ]; % should be faster than concatenating to third dimension
                        
        end
    end
    
    
    fprintf(repmat('\b',1,lineLength))
end

%% EMPIRICAL DIFFERENCE LOW
powLowDiff = cat(3,powLowDiff{:});
powLowDiff = squeeze(nanmean(powLowDiff,3));

%% HIGH
powHighDiff = cat(3,powHighDiff{:});
powHighDiff = squeeze(nanmean(powHighDiff,3));

 %% PERMUTATION TEST
 nperm = 100;
powLowPerm = zeros(nperm, 45, 1000);
for perm = 1:nperm
        lineLength = fprintf('%d permutations out of %d permutations done (%.2f%%).\n', perm, nperm, perm/nperm*100);

    powLowDiffComp = [];
    
    
    for comp = 1 : size(powLowComp,2)
        
        numIdx = powLowCompIdx(comp);
        oneComp = powLowComp{comp}(:,:,randperm(size(powLowComp{comp},3))); %% SHUFFLE
               
        
        powLowIdxPerm = squeeze(nanmean(oneComp(:,:,1:numIdx),3));
        powLowNdxPerm = squeeze(nanmean(oneComp(:,:,numIdx+1:end),3));
        
        
        powLowDiffComp = [powLowDiffComp, {powLowIdxPerm-powLowNdxPerm}];
        
    end
    
    powLowDiffComp = cat(3,powLowDiffComp{:});
    powLowPerm(perm,:,:) = squeeze(nanmean(powLowDiffComp,3));
end

save('powLowPerm', 'powLowPerm');

mean_h0 = squeeze(nanmean(powLowPerm,1));
std_h0  = squeeze(nanstd(powLowPerm,1));
zval = abs(norminv(0.05));

max_vals  = zeros(nperm,2);
bigIsland = zeros(nperm,1);
for perm = 1:nperm

    permMap = squeeze(powLowPerm(perm,:,:));
    permMap = (permMap-mean_h0)./std_h0;

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


% figure(5), clf
% histogram(max_vals(:),90)
% xlabel('Extremeiest H_0 mean difference values')
% ylabel('Count')


% find the lower threshold bound based on percentage
tmp = sort(max_vals(:,1));
thresh_lo = tmp( round(0.05*nperm) ); % SAME AS: prctile(tmp, 5)

% repeat for upper threshold bound
tmp = sort(max_vals(:,2));
thresh_hi = tmp( round((1-0.05)*nperm) );

figure(1); clf; hold on;
contmap = powLowDiff>=thresh_hi | powLowDiff <= thresh_lo;
imagesc(powLowDiff)

% and plot on top of the histogram
figure(5), hold on
plot([1 1]*thresh_lo,get(gca,'ylim'),'r--','linew',2)
plot([1 1]*thresh_hi,get(gca,'ylim'),'r--','linew',2)


% threshold real data
diffmap = squeeze(nanmean(powLowEmp,3));
zmap = squeeze(nanmean(powLowEmp,3));
zmap(zmap>thresh_lo & zmap<thresh_hi) = 0;

    %% OLD: POOLED OVER TRIALS INSTEAD OF LOOKING AT WITHIN (IDX-NDX) DIFFERENCE
    % allPowLowIdx = [];
    % allPowLowNdx = [];
    % allPowHiIdx  = [];
    % allPowHiNdx  = [];
    % for su = 1 : length(allSpks)
    %     lineLength = fprintf('%d spikes out of %d spikes done (%.2f%%).\n', su, length(allSpks), su/length(allSpks)*100);
    %
    %     idxTrlLw = allSpks(su).idxTrlSingLw;
    %     idxTrlHi = allSpks(su).idxTrlSingHi;
    %
    %     %% LOW
    %     if ~any(isnan(allSpks(su).idxTrlSingLw))
    %         allPowLowIdx = cat(3, allSpks(su).allPowLow(:,:, idxTrlLw), allPowLowIdx);
    %         allPowLowNdx = cat(3, allSpks(su).allPowLow(:,:,~idxTrlLw), allPowLowNdx);
    %     end
    %
    %     %% HIGH
    %     if ~any(isnan(allSpks(su).idxTrlSingHi))
    %         allPowHiIdx = cat(3, allSpks(su).allPowHigh(:,:, idxTrlHi), allPowHiIdx);
    %         allPowHiNdx = cat(3, allSpks(su).allPowHigh(:,:,~idxTrlHi), allPowHiNdx);
    %     end
    %
    %
    %     fprintf(repmat('\b',1,lineLength))
    % end
    %
    % %% DB NORMALIZE
    % % LOW IDX
    % for trl = 1:size(allPowLowIdx,3)
    %     for freq = 1:size(allPowLowIdx,1)
    %         tt = allPowLowIdx(freq,:,trl);
    %         allPowLowIdx(freq,:,trl) = 10*log10(tt./squeeze(nanmean(tt,2)));
    %     end
    % end
    %
    % % LOW NDX
    % for trl = 1:size(allPowLowNdx,3)
    %     for freq = 1:size(allPowLowNdx,1)
    %         tt = allPowLowNdx(freq,:,trl);
    %         allPowLowNdx(freq,:,trl) = 10*log10(tt./squeeze(nanmean(tt,2)));
    %     end
    % end
    %
    % % HIGH IDX
    % for trl = 1:size(allPowHiIdx,3)
    %     for freq = 1:size(allPowHiIdx,1)
    %         tt = allPowHiIdx(freq,:,trl);
    %         allPowHiIdx(freq,:,trl) = 10*log10(tt./squeeze(nanmean(tt,2)));
    %     end
    % end
    %
    % % HIGH NDX
    % for trl = 1:size(allPowHiNdx,3)
    %     for freq = 1:size(allPowHiNdx,1)
    %         tt = allPowHiNdx(freq,:,trl);
    %         allPowHiNdx(freq,:,trl) = 10*log10(tt./squeeze(nanmean(tt,2)));
    %     end
    % end
    % % save('precuepow_hiANDlo_raw.mat', 'allPowHiIdx', 'allPowHiNdx', 'allPowLowIdx', 'allPowLowNdx', '-v7.3')
    % % load('precuepow_hiANDlo_raw.mat')
    % diffLo = squeeze(nanmean(allPowLowIdx,3)) - squeeze(nanmean(allPowLowNdx,3));
    % diffHi = squeeze(nanmean(allPowHiIdx,3))  - squeeze(nanmean(allPowHiNdx,3));
    %
    % figure(4); clf;
    % imagesc(diffLo);
    % figure(5); clf;
    % imagesc(diffHi);
    
