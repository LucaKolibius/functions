%% THIS FUNCTION IS NOT USING MANUAL SPIKESORTING FOR THE SU
function datTrans4Simon_raw

addpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\functions');
% addpath(genpath('\\analyse4.psy.gla.ac.uk\Project0310\RDS\Common\mcode\tbx\chronux_2_11\spectral_analysis'))
addpath(genpath('\\analyse4.psy.gla.ac.uk\project0310\RDS\Common\mcode\tbx\chronux_2_11\spectral_analysis'));
addpath('\\analyse4.psy.gla.ac.uk\project0310\RDS\Common\mcode\utils\LFP'); % cleanLFPfromLineNoise
addpath('\\analyse4.psy.gla.ac.uk\project0310\RDS\Common\mcode\helper'); % CleanLineNoise
addpath(genpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\toolboxes\wave_clus-master'));


allSbj = dir('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\P*');
allSbj([1 2 4 5]) = [];

dsFS      = 1000;
dsTrlTime = linspace(-7,7,14001);
dsFs    = 1000;


for sbjNum = 1 : size(allSbj,1)
    cd([allSbj(sbjNum).folder, filesep, allSbj(sbjNum).name])
    allSesh = dir('S*');
    
    bidsID = sub_ID_conversion(allSbj(sbjNum).name,1);
    sbj    = sub_ID_conversion(bidsID,1);
        
    for seshNum = 1 : length(allSesh)
        
        sesh  = allSesh(seshNum).name;
        if strcmp(allSbj(sbjNum).name, 'P05') && strcmp(allSesh(seshNum).name, 'S5')
            continue % I couldn't split S5 from S4
        end
        
        if strcmp(allSbj(sbjNum).name, 'P13') && strcmp(allSesh(seshNum).name, 'S4')
            continue % cannot yet get the trigger to match up
        end
        
        %% CD TO SESSION FOLDER
        cd([allSesh(seshNum).folder, filesep, allSesh(seshNum).name]);
        seshDate = dir('2*'); cd(seshDate.name);
        seshDate = seshDate.name;
        
        %% FIND SAMPLING RATE
        SR = findSR;
        
        %% LOGFILE
        abc = dir('*LogFile_EMtask.txt'); if size(abc,1)>1; error('too many logfiles!'); end
        [trlENC, trlRET, rtENC, ~, missIdx, hitIdx] = loadLogsSimon([abc.folder, filesep, abc.name]);
        [~, hitsIdx, ~, ~, ~, retTrig, encTrig, ~, ~] = loadLogs([cd, filesep], 1); % GET ENC AND RET TRIGGER OF HITS
        
        encTrig = encTrig(:,1);
        retTrig = retTrig(:,1);
        
        retTrig         = sortrows(retTrig);
        allTrig         = [encTrig, ones(size(encTrig,1),1); retTrig, ones(size(retTrig,1),1)*2];
        
        RTs = rtENC;
        trlSel = [trlENC'; trlRET'];

        %% AUTOMAITC SPIKE SORTING
        % positive detection
        par = set_parameters();
        par.detection = 'pos';
        par.sr = SR;
        Get_spikes('all', 'parallel', true, 'par', par);
        movefile('*.mat', 'posDetectRaw'); % moves all .mat files that include the just detected spikes into the folder "posDetect"
        
        % % negative detection
        par.detection = 'neg';
        Get_spikes('all', 'parallel', true, 'par', par);
        movefile('*.mat', 'negDetectRaw'); % moves all .mat files that include the just detected spikes into the folder "negDetect"
        
        % positive clustering
        cd posDetectRaw
        Do_clustering('all', 'parallel', true, 'make_plots', false);
        
        % % negative clustering
        cd ..\negDetectRaw
        Do_clustering('all', 'parallel', true, 'make_plots', false);
        
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
            
        !!! CONTINUE HERE
        %% SPIKES
        for chan = 1 : length(chanLab)
            chanNam = chanLab{chan};
            
%             %% DETERMINE IF CURRENT WIRE HAS ANY SPIKES
%             if ~any(strcmp(chanNam,chanWithSpks)) % current channel does not have any spike on it
%                 sortedSpikesSEG{1,chan}=struct('newSpikeTimes',[],'assignedCluster',[],'wavf',[],...
%                     'num_clus',[],'SD',[],'SpikeTimesSeg',[], 'trl',[],'assignedClusterSeg',[]);
%                 wltCoeffs{1,chan} = [];
%                 continue
%             end
            
            %% GOES THROUGH ALL SPIKES IN THAT WIRE OF INTEREST
            allSpkTms = [];
            allCl     = [];
            wvf       = [];
            coeffs    = [];
            
                %% LOAD WAVESHAPE FROM DETECT FOLDER OR detect_newSR
                abc = dir(['times_*', chanNam, '.mat']);
                if size(abc,1)>1; error('Naming conflict before loading file'); end
                load(abc.name, 'spikes', 'cluster_class', 'par', 'inspk')
                                
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
            allCl = reSort(2,:);
            allSpkTms = reSort(1,:);
            
            sortedSpikesSEG{1,chan}=struct('newSpikeTimes',allSpkTms,'assignedCluster',allCl,'wavf',wvf,...
                'num_clus',length(sameWireSpks) ,'SD',[],'SpikeTimesSeg',allNewSpks(1,:), 'trl',allNewSpks(3,:),'assignedClusterSeg',allNewSpks(2,:));
            
            wltCoeffs{1,chan} = [];
        end % END OF LOOP OVER SU IN THAT SESSION
        
        %% SAVE LFP & SPK VARIABLES
        savepath = ['\\analyse4.psy.gla.ac.uk\Project0309\Luca\4simon_raw\', sbj, filesep, 'fVSpEM', filesep, seshDate, filesep];
        mkdir(savepath);
        save([savepath, sbj, '_fVSpEM_',seshDate,'_spkDataStimLockedSegmented.mat'], 'Fs', 'RTs', 'chanLab', 'hitIdx', 'missIdx', 'sortedSpikesSEG', 'trlENC', 'trlRET', 'trlSel', 'wltCoeffs')
        save([savepath, sbj, '_fVSpEM_',seshDate,'_lfpDataStimLockedSegmenteddownsampled.mat'], 'LFPseg', 'RTs', 'chanLab', 'dsFs', 'dsTrlTime', 'hitIdx', 'missIdx', 'trlENC', 'trlRET', 'trlSel')
    end
    cd(savepath) %% WILL PRODUCE AN ERROR FOR THE LAST SU IF IT IS SKIPPED
    
end % OF FUNCTION