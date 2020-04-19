%% crosscorrelation to rid some artifacts
clear all
cd Z:\Luca\data\P02\S2\2016-07-13_14-14-25\posDetect
allSpks = dir('times_CSC*'); % load all wires
tableTimestamps=table.empty;
subjID = 'P02_S2';

% pos spikes
% create a table with the spike timestamps for all wires, subdivided into
% cluster (includes 0 cluster with artifacts)
for ia = 1:size(allSpks,1)
    load(allSpks(ia).name, 'cluster_class'); % load one wire
    
    wirename = allSpks(ia).name;
    wirename(end-3:end)='';
    wirename(1:10)='';
    wirename(end+1:end+3)='Pos';
    
    allClust = max(cluster_class(:,1)); % how many cluster are in one wire
    
    % create table variable-names
    varname = cell.empty;
    for ib = 0:allClust
        varname(1,size(varname,2)+1) = cellstr([wirename, sprintf('CL%d', ib)]);
    end
    
    % create empty table with variable-names
    tempTable = cell2table(cell(1,size(varname,2)), 'VariableNames',varname);
    
    % insert spiketimes into table
    for ic = 0:allClust
        tempTable{1 , ic+1} = {cluster_class(cluster_class(:,1)==ic, 2)}; % ic+1 because it starts at "0"
    end
    
    tableTimestamps = [tableTimestamps tempTable];
    clear tempTable
end

%% negative spikes
cd Z:\Luca\data\P02\S2\2016-07-13_14-14-25\negDetect
allSpks = dir('times_CSC*'); % load all wires

% create a table with the spike timestamps for all wires, subdivided into
% cluster (includes 0 cluster with artifacts)
for ia = 1:size(allSpks,1)
    load(allSpks(ia).name, 'cluster_class'); % load one wire
    
    wirename = allSpks(ia).name;
    wirename(end-3:end)='';
    wirename(1:10)='';
    wirename(end+1:end+3)='Neg';
    
    allClust = max(cluster_class(:,1)); % how many cluster are in one wire
    
    % create table variable-names
    varname = cell.empty;
    for ib = 0:allClust
        varname(1,size(varname,2)+1) = cellstr([wirename, sprintf('CL%d', ib)]);
    end
    
    % create empty table with variable-names
    tempTable = cell2table(cell(1,size(varname,2)), 'VariableNames',varname);
    
    % insert spiketimes into table
    for ic = 0:allClust
        tempTable{1 , ic+1} = {cluster_class(cluster_class(:,1)==ic, 2)}; % ic+1 because it starts at "0"
    end
    
    tableTimestamps = [tableTimestamps tempTable];
    clear tempTable
end

%% find xc between all clusters
% make histogram so the delay between two cells is fixed
BW = 3.2; % 0.1ms bin width
startHist = 0; % start histogram
tempMax=[];
for ia = 1:size(tableTimestamps,2)
    tempMax(1,ia)=max(tableTimestamps{1,ia}{1});
end
endHist=[]; % end histogram
endHist=max(tempMax); % find last recorded spike

dt = [startHist: BW :endHist]; % 0.1ms bins / linear space for histogram
maxlag = 20; % 2ms

% preallocate deleteClust
deleteClust = zeros(7, size(tableTimestamps,2) *size(tableTimestamps,2));

for ia=1:size(tableTimestamps,2)
    a1=hist(tableTimestamps{1,ia}{1},dt);
    for ib=1:size(tableTimestamps,2)
        if ia==ib
            continue
        end
        
        % cross-correlation between two noise cluster are not interesting
        varname1 = tableTimestamps.Properties.VariableNames{ia};
        varname2 = tableTimestamps.Properties.VariableNames{ib};
        if strcmp(varname1(end-2:end),'CL0') && strcmp(varname2(end-2:end),'CL0')
            continue
        end
        
        disp([num2str(ia), ' - ', num2str(ib)]);
        a2=hist(tableTimestamps{1,ib}{1},dt);
        [xc,indx] = xcorr(a1, a2, maxlag, 'coeff'); % cross correlation
        
        % indexing clusters that exceed a specific correlation
        deleteClust(1,size(deleteClust,2)+1) = ia; % first cluster
        deleteClust(2,size(deleteClust,2))   = ib; % second cluster
        deleteClust(3,size(deleteClust,2))   = max(xc); % record highest cross-correlation
        temp = indx(xc==max(xc)); % record index where cross-correlation is maximal
        temp2 = temp(randi(size(temp,2))); % if I have two values for the peak position, choose one at random (will likely only happen for very small XC)
        deleteClust(4,size(deleteClust,2))   = temp2;
    end
end


% trim columns in which the first row is a 0 (comes from preallocation in L88
deleteClust(:,deleteClust(1,:)==0)=[];

% save variables
save(['deleteClust',subjID ,'.mat'], 'deleteClust');
save(['tableTimestamps', subjID, '.mat'], 'tableTimestamps');

% only work with the xc above 0.33
th = 0.33;
mindx = deleteClust(3,:)>=th;
deleteClust_th = deleteClust(:,mindx);

% save new deleteClust under a new variable name
save(['deleteClust_th',subjID,'.mat'],'deleteClust_th');

%% deal with reference induced spikes
counterSame = 0; % refers to same polarity
counterDiff = 0; % refers to different polarity
deleteClust_th(7,:) = 0;  % delete all row 7 tokens
finalDeleteClust = [];
for ix=1:size(deleteClust_th,2)
    
    % if this xc has been deleted in a previous step, continue
    if deleteClust_th(1,ix)==0
        continue
    end
    
    %not sure what this does
    if isempty(deleteClust_th(1,ix)) || isempty(deleteClust_th(2,ix))
        continue
    end
    
    if ix+6 > size(deleteClust_th,2) % otherwise I would index deleteClust_th outside its size and produce an error
        continue
    elseif ~deleteClust_th(1,ix)==deleteClust_th(1,ix+6) % if the same cluster does not have a significant xc with at least 6 other clusters, I do not need to investigate further
        continue
    end
    
    varname1 = tableTimestamps.Properties.VariableNames{deleteClust_th(1,ix)};
    varname2 = tableTimestamps.Properties.VariableNames{deleteClust_th(2,ix)};
    
    % if a new cluster is used as a xc base, reset counter
    if ix>1 % prevents error when ix=1
        if deleteClust_th(1,ix)~=deleteClust_th(1,ix-1)
            counterSame = 0;
            counterDiff = 0;
            deleteClust_th(7,:)=0;
        end
    end
    
    if strcmp(varname1(1:end-7), varname2(1:end-7)) % same bundle
        if ~strcmp(varname1(1:end-6), varname2(1:end-6)) % cannot be on the same wire
            if deleteClust_th(4,ix) == 0 % zero-lag
                if any(regexp(varname1, 'Neg')==regexp(varname2, 'Neg')) || any(regexp(varname1,'Pos')==regexp(varname2, 'Pos'))
                    counterSame = counterSame+1;
                    deleteClust_th(7,ix) = 2; % this is a temporary token indexing HP referencing for same polarity
                elseif any(regexp(varname1, 'Neg')==regexp(varname2,'Pos')) || any(regexp(varname1,'Pos')==regexp(varname2,'Neg'))
                    counterDiff = counterDiff+1;
                    deleteClust_th(7,ix) = 1; % this is a temporary token indexing HP referencing for polarity change
                end
            end
        end
    end
    
    % (1) HP wire spikes are reflected in one polarity in the original wire and 7 times in the other wires from the same bundle
    % in the opposite polarity. Here we only want to keep the original signal from the referencing wire
    if counterDiff==7 % 7 zero-lag xc to a different polarity. this xc is from the HP wire itself
        fprintf('%d is due to referencing against HP wire (case: HP wire NOT referenced against itself).',deleteClust_th(1,ix));
        
        finalDeleteClust = [finalDeleteClust, deleteClust_th(2,deleteClust_th(7,:)==1)];
        disp(sprintf('Deleting %.f cluster in loop #%.f', sum(deleteClust_th(7,:)==1), ix));
        
        % do not test the deleted cluster anymore
        % think more about this
        cluster_raus = deleteClust_th(2,deleteClust_th(7,:)==1); % cluster that will be deleted
        temp = ismember(deleteClust_th(1:2,:),cluster_raus); % does this cluster appear somewhere else as a first or second cluster?
        temp2 = temp(1,:)+temp(2,:);
        temp2=temp2/2;
        temp2=round(temp2);
        temp2=logical(temp2); % logical index which xc needs to go
        deleteClust_th(:,temp2)=0;
        
        deleteClust_th(7,:) = 0;  % delete all row 7 tokens
        
    elseif counterSame==6 % 6 zero-lag xc to the same polarity
        
        %                                         % (2) this might not be necessary anymore!
        %                                         if counterDiff==1 % also one zero-lag xc to a different polarity. the different polarity is the HP-wire signal
        %                                             disp([deleteClust_th(1,ix), 'is due to referencing against HP wire (case: HP wire NOT referenced against itself.'])
        %                                             unDelete(1,ix) = 1; % index cluster that should remain
        %                                             % delete all row 7 tokens
        
        % (3) same as (1), but the original HP is referenced against itself (and does not contain a signal anymore). Here
        % we only keep one of the 7 cluster in that bundle
        if counterDiff==0 % no zero-lag xc to a different polarity
            disp([deleteClust_th(1,ix), 'is possibly due to referencing against HP wire (case: HP wire referenced against itself).'])
            
            
            
            %% start visualisation
            % this needs to be visualized becasue it could be an artifact
            % occouring on the whole bundl
            beep
            flaggedCl = find(deleteClust_th(7,:)==2); % index of cluster that were flagged
            secondCl= deleteClust_th(2, flaggedCl); % number of 2nd cluster that were flagged
            firstCl = unique(deleteClust_th(1,flaggedCl)); % same for first
            
            varname1 = tableTimestamps.Properties.VariableNames{secondCl(1)};
            polarity = varname1(end-5:end-3);
            
            varname2 = tableTimestamps.Properties.VariableNames{secondCl(2)};
            varname3 = tableTimestamps.Properties.VariableNames{secondCl(3)};
            varname4 = tableTimestamps.Properties.VariableNames{secondCl(4)};
            varname5 = tableTimestamps.Properties.VariableNames{secondCl(5)};
            varname6 = tableTimestamps.Properties.VariableNames{secondCl(6)};
            keepCl = tableTimestamps.Properties.VariableNames{firstCl};
            
            cd ..
            if strcmp(polarity, 'Neg')
                cd negDetect
                disp('CD to negDetect');
            elseif strcmp(polarity, 'Pos')
                cd posDetect
                disp('CD to posDetect');
            end
            
            drawReferenceInduced(keepCl  , 2,4,1);
            drawReferenceInduced(varname1, 2,4,2);
            drawReferenceInduced(varname2, 2,4,3);
            drawReferenceInduced(varname3, 2,4,4);
            drawReferenceInduced(varname4, 2,4,5);
            drawReferenceInduced(varname5, 2,4,6);
            drawReferenceInduced(varname6, 2,4,7);
            % end visualisation
            
            % manual input whether the visualisation shows artefact or
            % reference induced
            temp = [];
            while isempty(temp)
                myResp = input('artefact or reference induced? (a/r) ', 's');
                temp = findstr(myResp, 'ar');
                if isempty(temp)
                    disp('Invalid response. Try again.');
                end
            end
            
            if strcmp(myResp,'a')
                deleteClust_th(7,:) = 0;  % delete all row 7 tokens
                continue % ? will be kicked out later, possibly with more (negative) clusters
            elseif strcmp(myResp,'r')
                disp(sprintf('Deleting %.f cluster in loop #%.f', sum(deleteClust_th(7,:)==2), ix));
                cluster_raus = deleteClust_th(2,deleteClust_th(7,:)==2); % cluster that will be deleted (only the second ones)
                finalDeleteClust = [finalDeleteClust, cluster_raus];
                temp = ismember(deleteClust_th(1:2,:),cluster_raus); % does this cluster appear somewhere else as a first or second cluster?
                temp2 = temp(1,:)+temp(2,:);
                temp2=temp2/2;
                temp2=round(temp2);
                temp2=logical(temp2); % logical index which xc needs to go
                deleteClust_th(:,temp2)=0;
            end
            
            deleteClust_th(7,:) = 0;  % delete all row 7 tokens
        end
    end
end

% importantly an artifact would create spikes in all wires at the same
% time, not in only a subset (as stated here)

% if I already
% discarded one of the clusters in one wire as an artifact in one
% wire i'll run into problems detecgting reference induced spikes

% hyper-polarization induced spikes
for ic=1:size(deleteClust_th,2)
    
    if isempty(deleteClust_th(1,ic)) || isempty(deleteClust_th(2,ic))
        continue
    end
    
    if deleteClust_th(1,ic)==0
        continue
    end
    
    disp(sprintf('XC %d of %d', ic, size(deleteClust_th,2)));
    varname1 = tableTimestamps.Properties.VariableNames{deleteClust_th(1,ic)};
    varname2 = tableTimestamps.Properties.VariableNames{deleteClust_th(2,ic)};
    
    % case when same microwire, but positive and negative (hyper-polarization
    % induced) == ONLY DELETE ONE
    if strcmp(varname1(1:end-6), varname2(1:end-6)) % same microwire
        if any(regexp(varname1, 'Neg')==regexp(varname2,'Pos')) || any(regexp(varname1,'Pos')==regexp(varname2,'Neg')) % polarity change
            
            % another approach would be to delete the cluster that comes after
            % (instead of looking at max amplitude)
            if ~deleteClust_th(4,ic) == 0 % hyper-polarization spikes cannot be lag-0
                
                % determine polarity of varname1 and cd to corresponding folder
                if regexp(varname1, 'Pos') % varname1 is negative
                    varname1pol = 'Pos';
                    cd posDetect
                elseif regexp(varname1, 'Neg') % varname1 is negative
                    varname1pol = 'Neg';
                    cd negDetect
                end
                
                % load spikedata of varname1
                load(['times_CSC_', varname1(1:end-6),'.mat'], 'spikes', 'cluster_class')
                clusterNum = str2double(varname1(end));
                indx2 = cluster_class(:,1)==clusterNum;
                
                % find amplitude of varname1
                clusterSpikes = spikes(indx2, :);
                clusterSpikes_avg = mean(clusterSpikes);
                clusterSpikes_max = max(abs(clusterSpikes_avg));
                clear spikes cluster_class
                
                % now we load varname2
                % switch to other folder
                if strcmp(varname1pol, 'Pos')
                    cd ..\negDetect
                elseif strcmp(varname1pol, 'Neg')
                    cd ..\posDetect
                end
                
                load(['times_CSC_', varname2(1:end-6),'.mat'], 'spikes', 'cluster_class')
                clusterNum = str2double(varname2(end));
                indx2 = cluster_class(:,1)==clusterNum;
                
                % find amplitude of varname2
                clusterSpikes = spikes(indx2, :);
                clusterSpikes_avg = mean(clusterSpikes);
                clusterSpikes_max2 = max(abs(clusterSpikes_avg));
                
                if clusterSpikes_max >= clusterSpikes_max2 % amplitude of varname1 > amplitude of varname2
                    if deleteClust_th(ic,4)<0 % max(xc) should be below 0 if time_x is before time_y
                        disp('Hyperpolarization induced spikes: Amplitude and XC are in accordance');
                        deleteClust_th(6,ic) = 2; % delete only second cluster / varname2
                    else
                        warning('Hyperpolarization induced spikes: Amplitude and XC are not in accordance. Breaking script.');
                        disp(ic);
                        beep
                        break
                    end
                else % amplitude of varname2 > amplitude of varname1
                    if deleteClust_th(ic,4)>0
                        disp('Hyperpolarization induced spikes: Amplitude and XC are in accordance');
                        deleteClust_th(6,ic) = 1; % delete only first cluster / varname1
                    else
                        warning('Hyperpolarization induced spikes: Amplitude and XC are not in accordance. Breaking script.');
                        disp(ic);
                        beep
                        break
                    end
                end
            end
        end
    end
end

%% case for artifacts occuring in multiple wires of the same bundle
% when within one bundle, various cluster on different MW ~xc with lag-0 it
% is likely an artifact
deleteClust_th(7,:) = 0; % reset row 7 marker; I use "3" as a marker in this part
goDelete = 0;
for ic=1:size(deleteClust_th,2)
    
    if deleteClust_th(1,ic)==0
        continue
    end
    
    disp(sprintf('XC %d of %d', ic, size(deleteClust_th,2)));
    varname1 = tableTimestamps.Properties.VariableNames{deleteClust_th(1,ic)};
    varname2 = tableTimestamps.Properties.VariableNames{deleteClust_th(2,ic)};
    
    counter3 = 0;
    if strcmp(varname1(1:end-7), varname2(1:end-7)) % both cluster are from the same bundle
        if deleteClust_th(4,ic) == 0 % lag-0
            counter3 = counter3+1;
            deleteClust_th(7,ic) = 3;
        end
    end
    
    % counter reset
    if ic > 1
        prevVarname1 = tableTimestamps.Properties.VariableNames{deleteClust_th(1,ic-1)};
        if ~strcmp(varname1(1:end-7),prevVarname1(1:end-7))
            counter3 = 0;
            goDelete = 1;
        end
    end
    
    % decision & marker reset
    % IMPORTANT: i compare it against negative AND positive polarities! So its
    % 15 comparisons. But they are not independent..
    if counter3 >= 3 && goDelete==1
        deleteClust_th(6,deleteClust_th(7,:)==3) = 5; % mark both cluster as to be deleted if they have a "3" mark in row 7
        deleteClust_th(7,:) = 0; % reset marker
        goDelete = 0;
    end
end
%% if any CL0 ~xc w/ another cluster => delete other cluster
for ic=1:size(deleteClust_th,2)
    
    if deleteClust_th(1,ic)==0
        continue
    end
    
    disp(sprintf('XC %d of %d', ic, size(deleteClust_th,2)));
    varname1 = tableTimestamps.Properties.VariableNames{deleteClust_th(1,ic)};
    varname2 = tableTimestamps.Properties.VariableNames{deleteClust_th(2,ic)};
    
    if strcmp(varname1(end),'0') || strcmp(varname2(end),'0')
        deleteClust_th(6,ix) = 5; % delete both
    end
end

% update finalDeleteClust
for iz = 1:size(deleteClust_th,2)
    if deleteClust_th(6,iz)==0
        continue
    elseif deleteClust_th(6,iz)==1
        finalDeleteClust = [finalDeleteClust, deleteClust_th(1,iz)];
    elseif deleteClust_th(6,iz)==2
        finalDeleteClust = [finalDeleteClust, deleteClust_th(2,iz)];
    elseif deleteClust_th(6,iz)==5
        finalDeleteClust = [finalDeleteClust, deleteClust_th(1,iz), deleteClust_th(2,iz)];
    end
end

save(['finalDeleteClust', subjID, '.mat'], 'finalDeleteClust');
tableTimestamps_postXC=tableTimestamps;
tableTimestamps_postXC(:, [finalDeleteClust]) = []; % delete clusters from table
save(['tableTimestamps_postXC_', subjID, '.mat'], 'tableTimestamps_postXC');


% %% use deleteClust to change respective cluster to CL0
% % do a backup first!
% % move all the times_CSC_<wirename>.mat into one folder
% positive
% cd Z:\Luca\data\P02\S4\2016-07-10_18-41-47\posDetect
% movefile('times_CSC_*.mat', 'Z:\Luca\data\P02\S4\2016-07-10_18-41-47\backup\posDetect'); % moves all .mat files that include the just detected spikes into the folder "posDetect"
%
% % work the magic
% for ie = 1:size(finalDeleteClust,2)
%             varname = tableTimestamps.Properties.VariableNames{finalDeleteClust(ie)};
%             wire = varname(1:end-6); % test this!
%             cluster = varname(end);
%             load(['times_CSC',wire,'.mat']) % test this
%             cluster_class(cluster_class(:,1)==cluster,2) = 0;
%             save(['times_CSC',wire,'.mat'], 'par', 'spikes', 'inspk', 'Temp', 'forced', 'gui_status', 'cluster_class', 'ipermut');
% end


% % negative
% cd Z:\Luca\data\P02\S4\2016-07-10_18-41-47\negDetect
% movefile('times_CSC_*.mat', 'Z:\Luca\data\P02\S4\2016-07-10_18-41-47\backup\negDetect'); % moves all .mat files that include the just detected spikes into the folder "posDetect"
% for ig = 1:size(finalDeleteClust,2)
%             varname = tableTimestamps.Properties.VariableNames{finalDeleteClust(ie)};
%             wire = varname(1:end-6); % test this!
%             cluster = varname(end);
%             load(['times_CSC',wire,'.mat']) % test this
%             cluster_class(cluster_class(:,1)==cluster,2) = 0;
%             save(['times_CSC',wire,'.mat'], 'par', 'spikes', 'inspk', 'Temp', 'forced', 'gui_status', 'cluster_class', 'ipermut');
% end

%% Visualisation of XC
% a1=hist(tableTimestamps{1,deleteClust(1,ic)}{1},dt);
% a2=hist(tableTimestamps{1,deleteClust(2,ic)}{1},dt);
% [xc,indx] = xcorr(a1, a2, maxlag, 'coeff');
% mHandle = plot(indx , xc);
% hold on
%
% % figure parameters
% title(['Cross-correlation between ', char(10), varname1,' and ', varname2]);
% xticks([indx(1):1:indx(end)]);
% xlim([indx(1) indx(end)]);
% xlabel(sprintf('Maximal overlap: %.3f', max(xc)));
% yticks([-0.1:0.1:1]);
% ylim([-0.05, 1]);
%
% % red line symbolizing the threshold / cut-off value of xc
% plot([indx(1),indx(end)], [0.45, 0.45], 'color', 'r');
% hold off

fireTh = endHist/320000; % frequency TH for cluster
clearvars -except fireTH tableTimestamps_postXC
p2d = cd;
p2d(end+1)='\';
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, hitsIdx, missIdx, allTrials, sessionNum, retTrigger, encTrigger]=loadLogs(p2d);

% spikenumber tables
encSpikeNumber_cueLocked=encSpiketimes_cueLocked;
encSpikeNumber_respLocked=encSpiketimes_cueLocked;
retSpikeNumber_cueLocked=encSpiketimes_cueLocked;
retSpikeNumber_respLocked=encSpiketimes_cueLocked;

%% timewindows
timeWindow=[];
timeWindow.encCueLocked=[-1 5];
timeWindow.encRespLocked=[-2 1];
timeWindow.retCueLocked=[-1 2];
timeWindow.retRespLocked=[-2 1];

clNames = {};
for iy = 1:size(tableTimestamps_postXC,2)
    varname = tableTimestamps_postXC.Properties.VariableNames{iy};
    
    if strcmp(varname(end), '0')
        disp(varname);
        continue
    end
    
    clNames(1, size(clNames,2)+1) = cellstr(varname);
    
    a = table2array(tableTimestamps_postXC{1,iy});
    a = a/32000; % samples to secs
    
    % add spiketimes of the cluster to the table
    % Encoding + CueLocked
    locking=1; % for CueLocked
    encSpiketimes_cueLocked{size(encSpiketimes_cueLocked,1)+1,:}=insertSpiketimes(encTrigger,a,locking,timeWindow.encCueLocked);
    
    % Encoding + RespLocked
    locking=3; % for RespLocked
    encSpiketimes_respLocked{size(encSpiketimes_respLocked,1)+1,:}=insertSpiketimes(encTrigger,a,locking,timeWindow.encRespLocked);
    
    % Retrieval + CueLocked
    locking=1;
    retSpiketimes_cueLocked{size(retSpiketimes_cueLocked,1)+1,:}=insertSpiketimes(retTrigger,a,locking,timeWindow.retCueLocked);
    
    % Retrieval + RespLocked
    locking=3;
    retSpiketimes_respLocked{size(retSpiketimes_respLocked,1)+1,:}=insertSpiketimes(retTrigger,a,locking,timeWindow.retRespLocked);
    
end

beep
artefacts_man = [12:20];
clNames(artefacts_man-2) = [];
save(['artefacts_man', subjID, '.mat'], 'artefacts_man');
save(['clNames', subjID, '.mat'], 'clNames');
encSpiketimes_cueLocked(artefacts_man,:) = [];
encSpiketimes_respLocked(artefacts_man,:) = [];
retSpiketimes_cueLocked(artefacts_man,:) = [];
retSpiketimes_respLocked(artefacts_man,:) = [];

save(['encSpiketimes_cueLocked', subjID, '.mat'],'encSpiketimes_cueLocked');
save(['encSpiketimes_respLocked', subjID, '.mat'],'encSpiketimes_respLocked');
save(['retSpiketimes_cueLocked', subjID, '.mat'], 'retSpiketimes_cueLocked');
save(['retSpiketimes_respLocked', subjID, '.mat'],'retSpiketimes_respLocked');

% transform spiketimes into number of spikes
% consider preallocating for speed
% no idea why i did num2cell here..
for ia=1:size(encSpiketimes_cueLocked,2)
    for ib=3:size(encSpiketimes_cueLocked,1)
        encSpikeNumber_cueLocked{ib,ia}=num2cell(cellfun(@length,encSpiketimes_cueLocked{ib,ia})/2); % because the time window is 2x as long
        encSpikeNumber_respLocked{ib,ia}=num2cell(cellfun(@length,encSpiketimes_respLocked{ib,ia}));
        retSpikeNumber_cueLocked{ib,ia}=num2cell(cellfun(@length,retSpiketimes_cueLocked{ib,ia}));
        retSpikeNumber_respLocked{ib,ia}=num2cell(cellfun(@length,retSpiketimes_respLocked{ib,ia}));
    end
end

% cueLocked Hits
[tempCell_enc, tempCell_ret]=mk_tempCell(allTrials, encSpikeNumber_cueLocked, retSpikeNumber_cueLocked, 'hit', 0);

clusterMean=[];
clusterStd=[];
for i=1:size(tempCell_enc,1)
    clusterMean(i,1)=mean([tempCell_enc{i,:} tempCell_ret{i,:}]);
    clusterStd(i,1)=std([tempCell_enc{i,:} tempCell_ret{i,:}]);
end

RSAhits_cueLocked=normSpikeNumber(tempCell_enc, tempCell_ret, clusterMean, clusterStd);

% cueLocked Misses
[tempCell_enc, tempCell_ret]=mk_tempCell(allTrials, encSpikeNumber_cueLocked, retSpikeNumber_cueLocked, 'miss', 0);
RSAmiss_cueLocked=normSpikeNumber(tempCell_enc, tempCell_ret, clusterMean, clusterStd);

% respLocked Hits
[tempCell_enc, tempCell_ret]=mk_tempCell(allTrials, encSpikeNumber_respLocked, retSpikeNumber_respLocked, 'hit', 0);
RSAhits_respLocked=normSpikeNumber(tempCell_enc, tempCell_ret, clusterMean, clusterStd);

% respLocked Misses
[tempCell_enc, tempCell_ret]=mk_tempCell(allTrials, encSpikeNumber_respLocked, retSpikeNumber_respLocked, 'miss', 0);
RSAmiss_respLocked=normSpikeNumber(tempCell_enc, tempCell_ret, clusterMean, clusterStd);

% plotting
plotRSA(RSAhits_respLocked, RSAmiss_respLocked, RSAhits_cueLocked, RSAmiss_cueLocked, hitsIdx, missIdx);

%% number of FF/PP/FP for hits & misses
% hit
numFF_hit=0;
for i=1:size(encSpikeNumber_cueLocked,2)
    if strcmp(encSpikeNumber_cueLocked{1,i}, 'ff') && strcmp(encSpikeNumber_cueLocked{2,i}, 'hit')
        numFF_hit=numFF_hit+1;
    end
end

numPP_hit=0;
for i=1:size(encSpikeNumber_cueLocked,2)
    if strcmp(encSpikeNumber_cueLocked{1,i}, 'pp') && strcmp(encSpikeNumber_cueLocked{2,i}, 'hit')
        numPP_hit=numPP_hit+1;
    end
end
numPP_hit=numPP_hit+numFF_hit;

numFP_hit=0;
for i=1:size(encSpikeNumber_cueLocked,2)
    if strcmp(encSpikeNumber_cueLocked{1,i}, 'fp') && strcmp(encSpikeNumber_cueLocked{2,i}, 'hit')
        numFP_hit=numFP_hit+1;
    end
end
numFP_hit=numFP_hit+numPP_hit;
numHit=[numFF_hit numPP_hit numFP_hit];

% miss
numFF_miss=0;
for i=1:size(encSpikeNumber_cueLocked,2)
    if strcmp(encSpikeNumber_cueLocked{1,i}, 'ff') && strcmp(encSpikeNumber_cueLocked{2,i}, 'miss')
        numFF_miss=numFF_miss+1;
    end
end

numPP_miss=0;
for i=1:size(encSpikeNumber_cueLocked,2)
    if strcmp(encSpikeNumber_cueLocked{1,i}, 'pp') && strcmp(encSpikeNumber_cueLocked{2,i}, 'miss')
        numPP_miss=numPP_miss+1;
    end
end
numPP_miss=numPP_miss+numFF_miss;

numFP_miss=0;
for i=1:size(encSpikeNumber_cueLocked,2)
    if strcmp(encSpikeNumber_cueLocked{1,i}, 'fp') && strcmp(encSpikeNumber_cueLocked{2,i}, 'miss')
        numFP_miss=numFP_miss+1;
    end
end
numFP_miss=numFP_miss+numPP_miss;
numMiss=[numFF_miss numPP_miss numFP_miss];

%% GLM
% analyse RSA
% between category is 3
RSA_mask=zeros(size(hitsIdx,1));
RSA_mask=RSA_mask+3;

% within category is 2
for x=1:size(hitsIdx,1)
    for y=1:size(hitsIdx,1)
        if x <=numHit(1) && y<=numHit(1)
            RSA_mask(y,x)=2;
        elseif x>numHit(1) && x<=numHit(2) && y>numHit(1) && y<=numHit(2)
            RSA_mask(y,x)=2;
        elseif x>numHit(2) && y>numHit(2)
            RSA_mask(y,x)=2;
        end
    end
end

% main diagonal is 1
for x=1:size(hitsIdx,1)
    for y=1:size(hitsIdx,1)
        if x==y
            RSA_mask(y,x)=1;
        end
    end
end

% GLM for hits + respLocked
glmY=[];
glmX=zeros(size(hitsIdx,1),2);
for i=1:size(hitsIdx,1)*size(hitsIdx,1)
    glmY(i) = RSAhits_respLocked(i);
    if RSA_mask(i)==1 % main diagnonal
        glmX(i,1)=1;
        glmX(i,2)=1;
    elseif RSA_mask(i)==2
        glmX(i,2)=1;
    end
end

% [RSAcoeff_hitsResp, RSAdev_hitsResp, RSAstats_hitsResp]=glmfit(glmX,glmY,'normal');
RSA_hitsResp =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)

% GLM for hits + cueLocked
glmY=[];
glmX=zeros(size(hitsIdx,1),2);
for i=1:size(hitsIdx,1)*size(hitsIdx,1)
    glmY(i)=RSAhits_cueLocked(i);
    if RSA_mask(i)==1 % main diagnonal
        glmX(i,1)=1;
        glmX(i,2)=1;
    elseif RSA_mask(i)==2
        glmX(i,2)=1;
    end
end

% [RSAcoeff_hitsCue, RSAdev_hitsCue, RSAstats_hitsCue] = glmfit(glmX,glmY,'normal');
RSA_hitsCue =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)

% Misses
% between category is 3
% i can probably just reuse the old RSA_mask variable
RSA_mask=zeros(size(missIdx,1));
RSA_mask=RSA_mask+3;

% within category is 2
for x=1:size(missIdx,1)
    for y=1:size(missIdx,1)
        if x <=numMiss(1) && y<=numMiss(1)
            RSA_mask(y,x)=2;
        elseif x>numMiss(1) && x<=numMiss(2) && y>numMiss(1) && y<=numMiss(2)
            RSA_mask(y,x)=2;
        elseif x>numMiss(2) && y>numMiss(2)
            RSA_mask(y,x)=2;
        end
    end
end

% main diagonal is 1
for x=1:size(missIdx,1)
    for y=1:size(missIdx,1)
        if x==y
            RSA_mask(y,x)=1;
        end
    end
end

glmY=[];
glmX=zeros(size(missIdx,1),2);
for i=1:size(missIdx,1)*size(missIdx,1)
    glmY(i)=RSAmiss_cueLocked(i);
    if RSA_mask(i)==1 % main diagnonal
        glmX(i,1)=1;
        glmX(i,2)=1;
    elseif RSA_mask(i)==2
        glmX(i,2)=1;
    end
end

% [RSAcoeff_missCue, RSAdev_missCue, RSAstats_missCue] = glmfit(glmX,glmY,'normal');
RSA_missCue =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)

glmY=[];
glmX=zeros(size(missIdx,1),2);
for i=1:size(missIdx,1)*size(missIdx,1)
    glmY(i)=RSAmiss_respLocked(i);
    if RSA_mask(i)==1 % main diagnonal
        glmX(i,1)=1;
        glmX(i,2)=1;
    elseif RSA_mask(i)==2
        glmX(i,2)=1;
    end
end

% [RSAcoeff_missResp, RSAdev_missResp, RSAstats_missResp] = glmfit(glmX,glmY,'normal');
RSA_missResp =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)


save(['glmfit', subjID, '.mat'], 'RSA_hitsCue', 'RSA_hitsResp', 'RSA_missCue', 'RSA_missResp');