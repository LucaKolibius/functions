%cd('\\analyse4\Project0310\iEEG_DATA\MICRO\P07\fVSpEM\2017-05-06_19-13-42')
%cd('\\analyse4.psy.gla.ac.uk\project0309\Luca\4simon2\P07\fVSpEM\2017-05-06_19-13-42')
%% Spike Daten segmentiert (-7 bis 7 Sekunden um Cue Onset).
%  Das sind Files mit der Extension “_lfpDataStimLockedSegmenteddownsampled.mat”.
%
% Darin sind dann folgende Variablen enthalten:
% sortedSpikesSEG: pro channel spikes, spike times, waveshape, etc. das ist der output von Waveclus, das solltest du in etwa kennen. Wichtig ist dass die spike times hier in Sekunden angegeben werden, relativ zum Cue onset, in einem Intervall zwischen -7 und 7 Sekunden.
 
addpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\functions');
% addpath(genpath('\\analyse4.psy.gla.ac.uk\Project0310\RDS\Common\mcode\tbx\chronux_2_11\spectral_analysis'))
addpath(genpath('\\analyse4.psy.gla.ac.uk\project0310\RDS\Common\mcode\tbx\chronux_2_11\spectral_analysis'));
addpath('\\analyse4.psy.gla.ac.uk\project0310\RDS\Common\mcode\utils\LFP'); % cleanLFPfromLineNoise
addpath('\\analyse4.psy.gla.ac.uk\project0310\RDS\Common\mcode\helper'); % CleanLineNoise

clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat', 'allSpks');
allBids  = {allSpks.bidsID};
allSesh  = {allSpks.sesh};
allWire  = {allSpks.wirename};

prevID = [];
dsFs    = 1000;
for su = [284:625 1:283] %1:length(allSpks)
    clearvars -except allBids allSesh allWire prevID dsFs su allSpks
    cd('\\analyse4.psy.gla.ac.uk\project0309\Luca\data')
    bidsID          = allSpks(su).bidsID;
    sesh            = allSpks(su).sesh;
    sr              = allSpks(su).initSR;
    
    %     % only hits
    %     encTrig         = allSpks(su).encTrigger(allSpks(su).hitsIdx,[2]); % in seconds
    %     retTrig         = allSpks(su).retTrigger(allSpks(su).hitsIdx,[2]); % in seconds
    
    encTrig         = allSpks(su).encTrigger(:,[1]); % in seconds
    retTrig         = allSpks(su).retTrigger(:,[1]); % in seconds
    
    retTrig         = sortrows(retTrig);
    allTrig         = [encTrig, ones(size(encTrig,1),1); retTrig, ones(size(retTrig,1),1)*2];
%     allTrig         = sortrows(allTrig,1); % at least for the LFP Fred first has all enc then all ret trials
    
    
    %     hitIdx          = allSpks(su).hitsIdx;
    %     missIdx         = ones(1,size(allSpks(su).encTrigger,1));
    %     missIdx(hitIdx) = 0;
    %     missIdx         = find(missIdx)';
    sameSesh        = strcmp(bidsID, allBids) & strcmp(sesh, allSesh);
    sameSesh        = find(sameSesh);
        
    %% SKIP REPEATS
    curID = [bidsID, regexprep(sesh, 'S1b', 'S1')];
    if strcmp(prevID, curID) % same patient + session
        continue
    else
        prevID = curID;
    end
    
    %% LFP
    load(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data\microLFP\with spk int\',bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_SPKINT.mat'], 'data')
    chanLab = data.label';
    
%     % GET OUT LINENOISE
%     step = 3;  % 3 iterations
%     s    = 15; % 15 second moving time window
%     for chan = 1:size(chanLab,2)
%         lfp = data.trial{1}(chan,:);
%         [lfpClean] = cleanLFPfromLineNoise(lfp,1000,step,s);
%         data.trial{1}(chan,:) = lfpClean;
%     end
    
    for chan = 1 : size(chanLab,2)
        lfp = data.trial{1}(chan,:); %7s before and after for each trial
        
        for trl = 1:size(allTrig,1)
            trg = round(allTrig(trl,1)*1000);
            LFPseg{1,chan}(:,trl) = lfp(trg-7000:trg+7000);
        end
    end
    
    dsFS      = 1000;
    dsTrlTime = linspace(-7,7,14001);
    
    %% LOGFILE
    sbj = sub_ID_conversion(bidsID,1);
    cd(sbj); cd(sesh); seshDate = dir('2*'); cd(seshDate.name);
    seshDate = seshDate.name;
    seshFold = cd;
    
    abc = dir('*LogFile_EMtask.txt'); if size(abc,1)>1; error('too many logfiles!'); end
    [trlENC, trlRET, rtENC, ~, missIdx, hitIdx] = loadLogsSimon([abc.folder, filesep, abc.name]);
    RTs = rtENC;
    trlSel = [trlENC'; trlRET'];
    
    %% SPKS
    chanWithSpks = {allSpks(sameSesh).wirename};
    for chan = 1 : length(chanLab)
        chanNam = chanLab{chan};
        
        %% DETERMINE IF CURRENT WIRE HAS ANY SPIKES
        if ~any(strcmp(chanNam,chanWithSpks)) % current channel does not have any spike on it
            sortedSpikesSEG{1,chan}=struct('newSpikeTimes',[],'assignedCluster',[],'wavf',[],...
                'num_clus',[],'SD',[],'SpikeTimesSeg',[], 'trl',[],'assignedClusterSeg',[]);
            wltCoeffs{1,chan} = [];
            continue
        end
        
        %% GOES THROUGH ALL SPIKES IN THAT WIRE OF INTEREST
        allSpkTms = [];
        allCl     = [];
        wvf       = [];
        coeffs    = [];
        
        sameWireSpks    = strcmp(bidsID, allBids) & strcmp(sesh, allSesh) & strcmp(chanNam, allWire); % finds all spks on a wire
        sameWireSpks    = find(sameWireSpks);

        for wtWire = 1 : length(sameWireSpks)                                  % loop over SU in that same session
            suWire          = sameWireSpks(wtWire);
            wirename        = allSpks(suWire).wirename;
            isPos           = allSpks(suWire).isPos;
            curSU           = allSpks(suWire).su;
            spks            = allSpks(suWire).spks';                           % in 1000 hz
            spks            = spks/1000;                                       % in seconds
            Fs              = allSpks(suWire).initSR;
            
            switch isPos
                case 0
                    posNeg = 'neg';
                case 1
                    posNeg = 'pos';
            end
            
            %% GET WAVESHAPE
            cd(seshFold);
            if exist('detect_newSR', 'file') == 7
                cd('detect_newSR');
            else
                cd([posNeg, 'Detect']); % go to posDetect or negDetect
            end
            
            %% LOAD WAVESHAPE FROM DETECT FOLDER OR detect_newSR
            abc = dir(['times_*', wirename, '.mat']);
            if size(abc,1)>1; error('Naming conflict before loading file'); end
            load(abc.name, 'spikes', 'cluster_class', 'par', 'inspk')
            
            clIdx        = cluster_class(:,1) == curSU;
            numSpksVar   = sum(clIdx);
            numSpks      = size(spks,2);
            
            %% IF SPIKETIMES ARE NOT CORRECT IT IS PROBABLY IN detect_newSR
            if ~isequal(par.sr, sr) | ~isequal(numSpksVar, numSpks)
                cd('..\detect_newSR2');
                abc = dir(['times_*', wirename,'_' posNeg, '.mat']);
                if size(abc,1)>1; error('Naming conflict before loading file'); end
                load(abc.name, 'spikes', 'cluster_class', 'par', 'inspk')
                
                % re-calculate number of spikes
                clIdx        = cluster_class(:,1) == curSU;
                numSpksVar   = sum(clIdx);
            end
            
            if ~isequal(par.sr, sr) | ~isequal(numSpksVar, numSpks)
                error('Still not load the correct waveclus file')
            end
            
            %% SAVING TO A GENERAL VARIABLE
            allSpkTms = [allSpkTms, spks];
            allCl     = [allCl, ones(1,numSpks)*curSU];
            wvf       = [wvf, spikes(clIdx,:)'];
%             coeffs    = [coeffs; inspk(clIdx,:)];
            
        end
        
        %% REORGANIZE INTO FRED'S STRUCTURE
        %% SORT ALLS SPIKETIMES WITHIN THAT WIRE BY TIMESTAMP
        reSort = [allSpkTms; allCl; wvf];
        reSort = reSort';
        reSort = sortrows(reSort,1);
        reSort = reSort';
        
        %% LOOP THROUGH TRIALS
        allNewSpks = [];
        for trl = 1:length(allTrig)
            idx = allTrig(trl,1)-7 <= reSort(1,:) & allTrig(trl,1)+7 >= reSort(1,:);
            
            newSpks         = [];
            newSpks(1,:)    = reSort(1,idx)-allTrig(trl,1); % segmented spiketimes
            newSpks(2,:)    = reSort(2,idx); % segmented cluster
            newSpks(3,:)    = trl*ones(1,sum(idx)); % segmented trl
            newSpks(4:67,:) = reSort(3:end,idx); % segmented wvf (not saved)
            
            allNewSpks = [allNewSpks, newSpks];
        end
        
        wvf = reSort(3:end,:)';
        sortedSpikesSEG{1,chan}=struct('newSpikeTimes',allSpkTms,'assignedCluster',allCl,'wavf',wvf,...
            'num_clus',length(sameWireSpks) ,'SD',[],'SpikeTimesSeg',allNewSpks(1,:), 'trl',allNewSpks(3,:),'assignedClusterSeg',allNewSpks(2,:));
        
        wltCoeffs{1,chan} = [];
    end % END OF LOOP OVER SU IN THAT SESSION
    
    %% SAVE LFP & SPK VARIABLES
    savepath = ['\\analyse4.psy.gla.ac.uk\Project0309\Luca\4simon2\', sbj, filesep, 'fVSpEM', filesep, seshDate, filesep];
    mkdir(savepath);
    save([savepath, sbj, '_fVSpEM_',seshDate,'_spkDataStimLockedSegmented.mat'], 'Fs', 'RTs', 'chanLab', 'hitIdx', 'missIdx', 'sortedSpikesSEG', 'trlENC', 'trlRET', 'trlSel', 'wltCoeffs')
    save([savepath, sbj, '_fVSpEM_',seshDate,'_lfpDataStimLockedSegmenteddownsampled.mat'], 'LFPseg', 'RTs', 'chanLab', 'dsFs', 'dsTrlTime', 'hitIdx', 'missIdx', 'trlENC', 'trlRET', 'trlSel')
end
cd(savepath)