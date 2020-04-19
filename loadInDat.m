% handles the load-in of data for visuSU
function [hitsIdx, subjID, spksEnc, spksRet, sensTrlsFFPP, allSpikes, animalCues, averageWS, su, ffTrls, ppTrls] = loadInDat(datapath, allSU, indexNeuron)
    cd(datapath)
    load(allSU(indexNeuron).name)
    
    % cd
    try
        cd X:/Luca/data
    catch
        cd /media/ldk898/rds-share/Luca/data
    end
    
    if ~exist('subjID') % if we only have the BIDS format, switch to old subject nomenclature
        subj = sub_ID_conversion(bidsID, 'yes');
        
        if strcmp(sesh, '1b')
            subj= 'P07ERL';
            sesh = 'S1b';
        end
        
        subjID = [subj, '_', sesh];
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
    abc = dir; cd(abc(3).name);
    
    % logfile
    p2d = cd;
    p2d(end+1)='\';
    sensTrlsFFPP = find(sensTrlsFFPP==1);
    [~, hitsIdx, missIdx, ~, ~, retTrigger, encTrigger, ~, animalCues] = loadLogs(p2d, sensTrlsFFPP);
    
    % load single units from processed variable
     % load single units from raw data (old)
    cd X:\Luca\data\allSbj
    loadVar = ['allSU_', bidsID,'_', sesh]; load(loadVar); % load tableTimestamps_postXC
    for ixx = 1:size(allSU,1)
        clusterSpikes = allSU{ixx,3}/32000;
        
        % adjusted for ENC(-1-113) RET(-1113)
        spksEnc(:, ixx) = insertSpiketimes2(encTrigger, clusterSpikes, [1 3], [-3 1])'; % 3 seconds prior to cue trigger until 1 second after response trigger
        spksRet(:, ixx) = insertSpiketimes2(retTrigger, clusterSpikes, [1 3], [-3 3])'; % 3 seconds prior to cue trigger until 3 seconds post response trigger
        
        allSpikes(ixx).waveshape = allSU{ixx,4};
        averageWS(ixx,:) = mean(allSpikes(ixx).waveshape,1);
    end
end
    
%     %% load single units from raw data (old)
%     cd advancedAnalysis\postXCrej
%     loadVar = dir('tableTimestamps_postXC_*'); load(loadVar.name); % load tableTimestamps_postXC
%     cd ../manualRej/
%     loadVar = dir('clNames*'); load(loadVar.name); % load clNames
%     cd ../..
%     clCount = 0;
%     allSpikes = [];
%     for i=1:size(tableTimestamps_postXC,2)
%         mVarname = tableTimestamps_postXC.Properties.VariableNames{1,i};
%         if ismember(mVarname, clNames)
%             clCount = clCount +1;
%             clusterSpikes = table2array(tableTimestamps_postXC{1,i})'; % load in all spiketimes from cluster mVarname
%             clusterSpikes = clusterSpikes/32000; % samples to sec
%             
%             % segment trials of this SU into trials
%             % old: here I inserted spikes in a time window surrounding one
%             % trigger
%             %             spksEnc(:, clCount) = insertSpiketimes(encTrigger, clusterSpikes, 1, [-2 5])'; % cue locked; 6s  %I could also generate a big file with all subj data in it (using a struct) or put enc and ret in one variable, but what for?
%             %             spksRet(:, clCount) = insertSpiketimes(retTrigger, clusterSpikes, 1, [-2 3])'; % cue locked; 4s
%             
%             % adjusted for ENC(0-213) RET(-1113)
%             spksEnc(:, clCount) = insertSpiketimes2(encTrigger, clusterSpikes, [1 3], [-2 0])'; % 2 seconds prior to cue trigger until response trigger
%             spksRet(:, clCount) = insertSpiketimes2(retTrigger, clusterSpikes, [1 3], [-3 3])'; % 3 seconds prior to cue trigger until 3 seconds post response trigger
%             
%             % get waveshape
%             demGate = cd; % where the ncs files are
%             
%             if ~isempty(regexp(mVarname, 'Pos', 'ONCE'))
%                 posNeg = 'posDetect';
%             elseif ~isempty(regexp(mVarname, 'Neg', 'ONCE'))
%                 posNeg = 'negDetect';
%             end
%             
%             cd(posNeg)
%             
%             temp1 = mVarname;
%             temp1(end-5:end) = [];
%             
%             % fixes a problem with the file names yet again.
%             if size(temp1,2) == 6
%                 temp1(1:2) = [];
%             end
%             temp = dir(['times_CSC_', temp1, '.mat']);
%             
%             % problem due to different file names. this fixes it. Its not
%             % pretty, but does the job.
%             if isempty(temp)
%                 temp = dir(['times_Micro_', temp1, '.mat']);
%             end
%         
%             load(temp.name, 'spikes', 'cluster_class');
%             
%             clNum = mVarname(end);
%             extrCl = cluster_class(:,1)==str2num(clNum);
%             
%             counterX = 0;
%             for ib = 1 : size(extrCl,1)
%                 if extrCl(ib,1) == 1
%                     counterX = counterX+1;
%                     allSpikes(clCount).waveshape(counterX, :) = spikes(ib,:);
%                 end
%             end
%             averageWS(clCount,:) = mean(allSpikes(clCount).waveshape,1);
%             cd(demGate)
%         end
%     end
% end