%% TUNING

clear; clc;
allSbj = dir('X:\Declan\tuning\sub-*');

% CYCLE THROUGH PARTICIPANTS
for sbj = 1 : length(allSbj)
    cd([allSbj(sbj).folder, filesep, allSbj(sbj).name])
    
    allSesh = dir([cd, filesep, '20*']);
    
    % CYCLE THROUGH SESSIONS
    for sesh = 1 : length(allSesh)
        try
            spks = [];
            cd([allSesh(sesh).folder, filesep, allSesh(sesh).name])
            
            %% DETERMINE SAMPLING FREQUENCY
            [FS] = getSamplFreq;
            
            cd('Spike Data\posDetect')
            posMW = dir('times_CSC_*.mat');
            cd('..\negDetect')
            allMW = [posMW; dir('times_CSC_*.mat')];
            
            %% Did Fred fuck up the MW labels in this session?
            allNams = {allMW.name};
            if any(cellfun(@(x) contains(x, '9'), allNams))
                shiteLabs = 1;
            else
                shiteLabs = 0;
            end
            
        catch
            sprintf('No "PosDetect" or "NegDetect" in Subj %d, Sesh %d', sbj, sesh);
            continue
        end
        
        
        % CYCLE THROUGH MW
        for mw = 1 : length(allMW)
            
            % LOAD IN SPIKETIMES
            load([allMW(mw).folder, filesep, allMW(mw).name], 'cluster_class');
            cluster_class(:,2) = round(cluster_class(:,2)/FS*1000);
            
            % NAME OF MICROWIRE
            % some microwires are - for whatever reason - pooled over a hemisphere (L1:L24, R1:24)
            if shiteLabs == 1
                lab = allMW(mw).name; % MW filename
                lab = lab(11:end-4); % MW identifier (e.g. L24)
                side = lab(1); % hemisphere (e.g. L)
                lab = str2num(lab(2:end)); % MW number (e.g. 24)
                lab = char(floor((lab-0.1)/8)+65); % create a new bundle identifier (1-8: A, 9:
                
                bundName = [side lab];
            else
                bundName = allMW(mw).name;
                bundName = bundName(11:end-5);
            end
           
            if  mw > size(posMW,1)
                posNeg = 'neg';
            else
                posNeg = 'pos';
            end
            
            % CYCLE THROUGH CLUSTER
            for cl = 1 : max(cluster_class(:,1))
                spikename = [allMW(mw).name(11:end-4), posNeg, num2str(cl)];
                spks = [spks; {cluster_class(cluster_class(:,1) == cl, 2)}, bundName, spikename];
            end
        end
        
        %% LOAD LFP DATA
        %         load('X:\Luca\sub_9_sess_4_RAW_INT_DS.mat', 'LFP_data');
        try
            cd('..\..\LFP Data')
            lfpFile = dir('sub_*_RAW_INT_DS.mat');
            load(lfpFile.name, 'LFP_data');
        catch
            sprintf('No LFP folder in Subj %d, Sesh %d', sbj, sesh);
            continue
        end
        %         cd('Artefact Rejection');
        %         artFile = dir('sub*_th-8.mat');
        %         load(artFile.name, 'del_sampleinfo');
        
        %% TRANSFORM SPIKE TIMES TO TIME SERIES
        %  GAUSSIAN KERNEL
        mlength   = [-2000:1:2000];
        mSigma    = 500;
        mKernel   = normpdf(mlength,0,mSigma);
        mKernel   = mKernel/max(mKernel); % normalize peak to 1
        plot(mlength, mKernel)
        
        % move from spiketimes into histogram bins
        dt        = LFP_data.time{1}-0.0005;                                         % shift to the left
        dt(1)     = 0;                                                               % don't go out of bounds
        dt(end+1) = dt(end)+0.001;                                                   % add last bin
        dt        = dt * 1000;                                                       % make into samples (this is pretty much a vectro 0:1:lenRec now);
        histTimes = cellfun(@(x) histcounts(x, dt), spks(:,1), 'UniformOutput', false);
        
        % convert these bincounts into a time series via gaussian convolution
        trialSer  = cellfun(@(x)conv(mKernel, x), histTimes, 'UniformOutput', false); % every row is a SU / every column is a trial
        trialSer  = vertcat(trialSer{:});
        
        % get rid of edges introduced by the convolution
        prune     = (length(mlength)-1) /2;
        trialSer  = trialSer(:,prune+1:end-prune);
        
        %% GET TRIAL INFO
        cd('..')
        [img, trigger] = loadLogs_tuningBT([cd, filesep]);
        
        %% SORT INTO BUNDLES
        allBund = unique(spks(:,2));
        
        % NEW: FOR CONCEPT CELLS
        % load in tCells.mat and image_cell
        try
            cd('SU Analysis');
            load('tCells.mat');
            load('allSUdata.mat', 'image_cell');
        catch
            sprintf('No "SU ANALYSIS" in Subj %d, Sesh %d', sbj, sesh);
            continue
        end
        
        for bund = 1 : length(allBund)
            
            % INDEX WHICH LFPs ARE RELEVANT FOR THE CURRENT STRUCT
            curBund = allBund(bund);
            LFPbund = cellfun(@(x) x(1:end-1), LFP_data.label, 'un', 0);
            bundIdx = strcmp(curBund, LFPbund);
            
            % GET ALL THE SPIKES FROM THE CURRENT BUNDLE
            spkBundIdx = strcmp(curBund, spks(:,2));
            
            nwDat              = [];
            nwDat.trial        = trialSer(spkBundIdx,:);
            nwDat.carrier      = LFP_data.trial{1}(bundIdx,:);
            nwDat.fsample      = 1000;
            nwDat.time         = LFP_data.time;
            nwDat.labelCarrier = LFP_data.label(bundIdx);
            nwDat.label        = spks(spkBundIdx,2);
            
            % OLD. get the trigger
            % nwDat.placesTrl    = trigger(strcmp(img, 'p'));
            % nwDat.facesTrl     = trigger(strcmp(img, 'f'));
            
            % some tuning indices are double, some are cells. here I make them all
            % cells so I can work with them
            for corCel = 1 : size(tCells,2)
                if isa(tCells(5,corCel), 'double')
                    tCells(t,corCel) = {tCells(t,corCel)};
                end
            end
            ccBund  = cellfun(@(x) x(1:end-5), tCells(1,:), 'un', 0); % the bundles on which all concept cells are
            ccIdx   = strcmp(curBund, ccBund);                        % index of all concept cells that are on the current bundle
            ccIdx   = find(ccIdx == 1);                               % find-index
            tunedIm = vertcat(tCells{5,ccIdx});                       %
            
            % save to variable
            nwDat.trigger = trigger;
            nwDat.tunedIm = tunedIm;
            nwDat.allIm = image_cell(:,2);
            
            save(['X:\Luca\code4cont\spkTS\', allSbj(sbj).name, '_S', num2str(sesh), '_', curBund{1}], 'nwDat');
            
        end  % END OF BUNDLE LOOP
        
    end % END OF SESSION LOOP
    
end % END OF SUBJECT LOOP