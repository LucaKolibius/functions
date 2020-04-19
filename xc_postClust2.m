function quickfix = xc_postClust2(subjID)

% P04S1 a positive spike correlated with the whole cluster of negatives
% (reference induced). The algorithm automatically selected all negative
% ones as artefacts, but then asked for manual input in a later stage
% because all negative artefacts correlated enough within each other to ask
% for "artefact or reference". Maybe saying "reference" here does not even
% do any damage (but it causes an error), because the relevant cluster are
% within "finalDelete" anyway?

%for the server
try
    cd X:/Luca/data
catch
    cd /media/ldk898/rds-share/Luca/data
end

mSubject = subjID(1:end-3);
mSession = subjID(end-1:end);
cd(mSubject)
cd(mSession)
abc = dir;
cd(abc(3).name)
cd advancedAnalysis
cd allXC
mCD = cd;
abc1 = dir('deleteClust_th*');
abc2 = dir('tableTimestamps*');
load(abc1.name)
load(abc2.name)

% just the names of the numbers in deleteClust_th
allNam={};
for ih = 1 : size(deleteClust_th,2)
    allNam{1,ih} = tableTimestamps.Properties.VariableNames{deleteClust_th(1,ih)};
    allNam{2,ih} = tableTimestamps.Properties.VariableNames{deleteClust_th(2,ih)};
end
save(['allNam', subjID, '.mat'], 'allNam');

%% deal with reference induced spikes
counterSame = 0; % refers to same polarity
counterDiff = 0; % refers to different polarity
deleteClust_th(7,:) = 0;  % delete all row 7 tokens
finalDeleteClust = [];
saveOC = [];
for ix = 1 : size(deleteClust_th,2)
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
        
        % save the OC from later deletion during multiple wires, same bundle AR
        saveOC = [saveOC, deleteClust_th(1,deleteClust_th(7,:)==1)];
        
        deleteClust_th(7,:) = 0;  % delete all row 7 tokens
        
    elseif counterDiff>7
        % This case: When I have reference induced spikes from a positive (neg)
        % wire and it is NOT referenced against itself, but a reference induced cluster
        % is split in two cluster
        % I guess this could be the OC... but in reality this will be CL2
        % and therefore a smaller cluster and rarely the OC / even a proper
        % cluster / likely a reference induced cluster split in two
        figure(2)
        varnameX = tableTimestamps.Properties.VariableNames{deleteClust_th(2,ix)};
        drawReferenceInduced(varnameX  , 2, 2, counterDiff-7);
        
        temp = [];
        while isempty(temp)
            myResp = input('Also delete this cluster? (y/n) ', 's');
            temp = findstr(myResp, 'yn');
            if isempty(temp)
                disp('Invalid response. Try again.');
            end
        end
        
        % it is not necessary to delete this cluster from deleteClust_th because
        % it should not be possible anymore that the OC gets flagged up during the
        % later wire-wide artefact rejection
        if strcmp(myResp,'y')
            finalDeleteClust = [finalDeleteClust, deleteClust_th(2,ix)];
        end
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
            disp([num2str(deleteClust_th(1,ix)), 'is possibly due to referencing against HP wire (case: HP wire referenced against itself).'])
            
            
            
            %% start visualisation
            % this needs to be visualized becasue it could be an artifact
            % occouring on the whole bundle
            % you could switch the counter to a sum(2-flags)
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
            keepCl   = tableTimestamps.Properties.VariableNames{firstCl};
            
            cd(mCD)
            if strcmp(polarity, 'Neg')
                cd ../../negDetect
                disp('CD to negDetect');
            elseif strcmp(polarity, 'Pos')
                cd ../../posDetect
                disp('CD to posDetect');
            end
            
            disp(['Visualizing: ', subjID]);
            close(gcf)
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
            beep
            temp = [];
            while isempty(temp)
                myResp = input('artefact or reference induced? (a/r) ', 's');
                temp = findstr(myResp, 'ar');
                if isempty(temp)
                    disp('Invalid response. Try again.');
                end
            end
            
            % if reference induced, which one to keep?
            if strcmp(myResp,'a')
                deleteClust_th(7,:) = 0;  % delete all row 7 tokens
                continue % ? will be kicked out later, possibly with more (negative) clusters
            elseif strcmp(myResp, 'r')
                temp = [];
                while isempty(temp)
                    myResp2 = input('Which cluster is the OR (origninal reference) that you keep? (1-8) ', 's');
                    temp = strfind('12345678', myResp2);
                    if isempty(temp)
                        disp('Invalid response. Try again.');
                    end
                end
            end
            
            if strcmp(myResp,'r')
                if strcmp(myResp2, '1') % myResp2 is 1
                    cluster_raus = secondCl; % kick out the other cluster that are reference induced
                    saveOC = [saveOC, firstCl]; % the first cluster is OC, save it
                elseif ~isempty(findstr(myResp2, '234567')) % myResp2 is 2:7
                    myResp2      = str2num(myResp2);
                    secondClm1   = secondCl;
                    secondClm1(myResp2) = [];
                    cluster_raus = [firstCl, secondClm1];
                    saveOC = [saveOC, secondCl(myResp2-1)];
                end
            end
            finalDeleteClust = [finalDeleteClust, cluster_raus];
            
            deleteClust_th(7,:) = 0;  % delete all row 7 tokens
        end
        
    elseif counterSame>6
        % This case: When I have reference induced spikes from a positive (neg)
        % wire and it is referenced against itself, but a reference induced cluster
        % is split in two cluster
        figure(2)
        varnameX = tableTimestamps.Properties.VariableNames{deleteClust_th(2,ix)};
        drawReferenceInduced(varnameX  , 2, 2, counterSame-6);
        
        temp = [];
        while isempty(temp)
            myResp = input('Also delete this cluster? (y/n) ', 's');
            temp = findstr(myResp, 'yn');
            if isempty(temp)
                disp('Invalid response. Try again.');
            end
        end
        
        % it is not necessary to delete this cluster from deleteClust_th because
        % it should not be possible anymore that the OC gets flagged up during the
        % later wire-wide artefact rejection
        if strcmp(myResp,'y')
            finalDeleteClust = [finalDeleteClust, deleteClust_th(2,ix)];
        end
        deleteClust_th(7,:) = 0;  % delete all row 7 tokens
    end
end
cd(mCD)

% importantly an artifact would create spikes in all wires at the same
% time, not in only a subset (as stated here)

%% hyper-polarization induced spikes
quickfix = [];
for ic=1:size(deleteClust_th,2)
    % these next three lines are debatable
    % we might have a proper spike in one cluster that then correlates with
    % a zero cluster of the opposite spike because
    
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
%             if ~deleteClust_th(4,ic) == 0 % hyper-polarization spikes cannot be lag-0
                
                % if I have a proper positive spike with a strong
                % hyperpolarization, I will pick it up in the negative
                % clustering and delete it, so it will be in cl0
                if strcmp(varname1(end),'0') || strcmp(varname2(end),'0')
                    quickfix = [quickfix, ic];
                    continue
                end
                
                % determine polarity of varname1 and cd to corresponding folder
                if regexp(varname1, 'Pos') % varname1 is negative
                    varname1pol = 'Pos';
                    cd ../../posDetect
                elseif regexp(varname1, 'Neg') % varname1 is negative
                    varname1pol = 'Neg';
                    cd ../../negDetect
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
                    cd ../negDetect
                elseif strcmp(varname1pol, 'Neg')
                    cd ../posDetect
                end
                
                load(['times_CSC_', varname2(1:end-6),'.mat'], 'spikes', 'cluster_class')
                clusterNum = str2double(varname2(end));
                indx2 = cluster_class(:,1)==clusterNum;
                
                % find amplitude of varname2
                clusterSpikes = spikes(indx2, :);
                clusterSpikes_avg = mean(clusterSpikes);
                clusterSpikes_max2 = max(abs(clusterSpikes_avg));
                
                if clusterSpikes_max >= clusterSpikes_max2 % amplitude of varname1 > amplitude of varname2
                    if deleteClust_th(4,ic)<0 % max(xc) should be below 0 if time_x is before time_y
                        disp('Hyperpolarization induced spikes: Amplitude and XC are in accordance');
                        deleteClust_th(6,ic) = 2; % delete only second cluster / varname2
                    else
                        warning('Hyperpolarization induced spikes: Amplitude and XC are not in accordance. Breaking script.');
                        disp(ic);
                        beep
                        break
                    end
                else % amplitude of varname2 > amplitude of varname1
                    if deleteClust_th(4,ic)>0
                        disp('Hyperpolarization induced spikes: Amplitude and XC are in accordance');
                        deleteClust_th(6,ic) = 1; % delete only first cluster / varname1
                    else
                        warning('Hyperpolarization induced spikes: Amplitude and XC are not in accordance. Breaking script.');
                        disp(ic);
                        beep
                        break
                    end
                end
%             end
        end
    end
end

cd(mCD);
%% case for artifacts occuring in multiple wires of the same bundle
% when within one bundle, various cluster on different MW ~xc with lag-0 it
% is likely an artifact
deleteClust_th(7,:) = 0; % reset row 7 marker; I use "3" as a marker in this part
counter3 = 0;
counter3th = 6;     % IMPORTANT: i compare it against negative AND positive polarities! So its 15 comparisons (but not independent)

for ic=1:size(deleteClust_th,2)
    
    if deleteClust_th(1,ic)==0
        continue
    end
    
    disp(sprintf('XC %d of %d', ic, size(deleteClust_th,2)));
    varname1 = tableTimestamps.Properties.VariableNames{deleteClust_th(1,ic)};
    varname2 = tableTimestamps.Properties.VariableNames{deleteClust_th(2,ic)};
    
    % counter reset
    if ic > 1
        % do this in numbers in deleteClust_th ?
        prevVarname1 = tableTimestamps.Properties.VariableNames{deleteClust_th(1,ic-1)};
        if ~strcmp(varname1(1:end-6),prevVarname1(1:end-6)) || ic == size(deleteClust_th,2)% change in wire or last comparison
            if counter3 >= counter3th
                disp('Artifacts on multiple wires.');
                deleteClust_th(6,deleteClust_th(7,:)==3) = 5;
            end
            deleteClust_th(7,:) = 0; % reset "3" marker
            counter3 = 0;
        end
    end
    
    if strcmp(varname1(1:end-7), varname2(1:end-7)) % both cluster are from the same bundle
        if ~strcmp(varname1(1:end-6), varname2(1:end-6)) % but not the same wire
            if deleteClust_th(4,ic) == 0 % lag-0
                counter3 = counter3+1;
                deleteClust_th(7,ic) = 3;
            end
        end
    end
    
end

%% if any CL0 ~xc w/ another cluster => delete other cluster
for ic=1:size(deleteClust_th,2)
    if ismember(ic , quickfix) % these XC are because I delete (say move to Cl0) spikes that are clearly positive spikes during the clustering of negative detected spikes
        continue
    end
    
    if deleteClust_th(1,ic)==0
        continue
    end
    
    disp(sprintf('XC %d of %d', ic, size(deleteClust_th,2)));
    varname1 = tableTimestamps.Properties.VariableNames{deleteClust_th(1,ic)};
    varname2 = tableTimestamps.Properties.VariableNames{deleteClust_th(2,ic)};
    
    a(ic)=size(deleteClust_th,2);
    
    if strcmp(varname1(end),'0') || strcmp(varname2(end),'0')
        deleteClust_th(6,ic) = 5; % delete both
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

mtemp = ismember(finalDeleteClust,saveOC); % OC is the original cluster that created reference-induced cluster
finalDeleteClust(mtemp) = []; % delete it from finalDelteClust

cd ..
mkdir('postXCrej');
cd postXCrej
save(['finalDeleteClust', subjID, '.mat'], 'finalDeleteClust');
tableTimestamps_postXC = tableTimestamps;
tableTimestamps_postXC(:, [finalDeleteClust]) = []; % delete clusters from table
save(['tableTimestamps_postXC_', subjID, '.mat'], 'tableTimestamps_postXC');

end