% spikesorting und spikedetection
% add all toolboxes
clear all
addpath(genpath('Z:\Common\fieldtrip-20170618')); % fieldtrip
ft_defaults;
addpath(genpath('Z:\Luca\functions\wave_clus-master')); % waveclus 3.0
addpath('Z:\Luca\functions'); % my functions
addpath('Z:\Luca\functions\Neuralynx_19012019'); % Neuralynx (the commons folder function doesnt work on my PC)
addpath(genpath('Z:\Luca\TREBER\Scripts'))

% cd 'Z:\Luca\data'
% [p2d] = 'Z:\Luca\data';
p2d = cd;
p2d(end+1)='\';

[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, hitsIdx, missIdx, allTrials, sessionNum, retTrigger, encTrigger]=loadLogs(p2d);

% spikenumber tables
encSpikeNumber_cueLocked=encSpiketimes_cueLocked;
encSpikeNumber_respLocked=encSpiketimes_cueLocked;
retSpikeNumber_cueLocked=encSpiketimes_cueLocked;
retSpikeNumber_respLocked=encSpiketimes_cueLocked;

% % pretrial baseline
% encSpikeNumbers_pretrial=double.empty;
% retSpikeNumbers_pretrial=double.empty;

%% raw signal
CSCfiles = dir([p2d,'*.ncs']);
par = set_parameters_Bham(32000); %% [in the original script instead of 32000 the brackets were empty]
clusterCount=0;
clusterCategory='';

for it = 1%:length(CSCfiles)
    disp(sprintf('Microwire: %d von %d',it,size(CSCfiles,1)));
%     fprintf('Microwire: %d von %d',it,size(CSCfiles,1));

    % I RESET dum1 AND dum2 HERE BECAUSE OF MEMORY CONSTRAINTS!!
    dum1 = cell(1,length(CSCfiles));
    dum2 = cell(1,length(CSCfiles));
    
    % define the FieldSelection variable and set the ExtracMode.
    FieldSelection(1) = 1;%timestamps
    FieldSelection(2) = 0;
    FieldSelection(3) = 0;%sample freq
    FieldSelection(4) = 0;
    FieldSelection(5) = 1;%samples
    ExtractHeader = 1;
    ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.
    
    % extract the channel name
    chanLab = CSCfiles(it).name;
    chanLab(regexp(chanLab,'_')) = [];
    chanLab(regexp(chanLab,'CSC'):regexp(chanLab,'CSC')+2) = [];
    chanLab(regexp(chanLab,'.ncs'):end) = [];
    
    % import raw data
    [timestampsCSC, dataSamplesCSC, hdrCSC] = Nlx2MatCSC(CSCfiles(it).name, FieldSelection, ExtractHeader, ExtractMode, []);
    
    % extract the scale factor
    chck = regexp(hdrCSC,'ADBitVolts');
    selIdx = [];
    for jt = 1:length(chck)
        selIdx(jt) = ~isempty(chck{jt});
    end
    selIdx = find(selIdx~=0);
    scalef = str2double(hdrCSC{selIdx}(min(regexp(hdrCSC{selIdx},'\d')):end));
    
    % flatten
    dataSamplesCSC = double(dataSamplesCSC(:))'; % raw data in samples
    dataSamplesCSC = dataSamplesCSC.*scalef.*1e6;
    
    % extract the sampling-frequency
    chck = regexp(hdrCSC,'SamplingFrequency');
    selIdx = [];
    for jt = 1:length(chck)
        selIdx(jt) = ~isempty(chck{jt});
    end
    selIdx = find(selIdx~=0);
    Fs = str2double(hdrCSC{selIdx}(min(regexp(hdrCSC{selIdx},'\d')):end));
    [TimeStampPerSample] = 1/Fs*1e6; % Fs: sampling rate => returns the timestamps recorded at each data sample
    
    % create a time vector of theCSCtime
    [CSCTime] = 0:1/Fs:(length( dataSamplesCSC )-1)/Fs; % in units of sec %creates 31.25ms bins
    
    %[SPKTime] = linspace(timestampsCSC(1),timestampsCSC(end) + 512*TimeStampPerSample,length(dataSamplesCSC)); % in units of micro-sec
    [SPKTime] = uint64(timestampsCSC(1):TimeStampPerSample:timestampsCSC(1)+(length(dataSamplesCSC)-1)*TimeStampPerSample);
    
    % spike detection
%     [~,spikeWaveforms,thr1,spikeTimestamps,noise_std_detect,noise_std_sorted]
%     = amp_detect(dataSamplesCSC,par); % old
    [spikeWaveforms,thr1,spikeTimestamps] = amp_detect(dataSamplesCSC,par);
    % Get_spikes does not even require loading with Nlx

    durWaveform = size(spikeWaveforms,2); % [was spikeWaveforms1,2] & = 64 samples
    wft = linspace(0,(durWaveform-1)/Fs,durWaveform); % time vector of waveform % linear spaced vector from a to b in c points / from 0 to 0.002 (spikewave length) in steps of 31.25 (sample freq)
    
    % TTL rejection (spikeWaveforms & spikestamps!)
    spikeTimestampsSec=spikeTimestamps./32000;
    
    % encCue
    timeWindow=[0 0.005];
    TTLidx=[];
    for ia=1:size(retTrigger,1)
        trialTTLIdx=find(spikeTimestampsSec>=encTrigger(ia,1)+timeWindow(1) & spikeTimestampsSec<=encTrigger(ia,1)+timeWindow(2));
        for ib=1:size(trialTTLIdx,2)
            TTLidx(size(TTLidx,2)+1)=trialTTLIdx(ib);
        end
    end
    
    % retCue
    timeWindow=[0 0.02];
    for ia=1:size(retTrigger,1)
        trialTTLIdx=find(spikeTimestampsSec>=retTrigger(ia,1)+timeWindow(1) & spikeTimestampsSec<=retTrigger(ia,1)+timeWindow(2));
        for ib=1:size(trialTTLIdx,2)
            TTLidx(size(TTLidx,2)+1)=trialTTLIdx(ib);
        end
    end
    
    % ret stim
    timeWindow=[-0.0015 0.02];
    for ia=1:size(retTrigger,1)
        trialTTLIdx=find(spikeTimestampsSec>=retTrigger(ia,2)+timeWindow(1) & spikeTimestampsSec<=retTrigger(ia,2)+timeWindow(2));
        for ib=1:size(trialTTLIdx,2)
            TTLidx(size(TTLidx,2)+1)=trialTTLIdx(ib);
        end
    end
    
    disp(sprintf('%d TTL spikes were removed from a total of %d spikes in this MW (%.2f%%)', size(TTLidx,2), size(spikeTimestampsSec,2), size(TTLidx,2)/size(spikeTimestampsSec,2)*100));
    spikeTimestampsSec=[];
    spikeTimestamps(TTLidx)=[];
    spikeWaveforms(TTLidx,:)=[];
    
    % SpikeSorting
    waveclus.spikes = spikeWaveforms;
    waveclus.index = spikeTimestamps;
    [dim] = size(spikeWaveforms);
    sortedSpikes.newSpikeTimes = [];
    sortedSpikes.assignedCluster = [];
    sortedSpikes.wavf = [];
    sortedSpikes.num_clus = [];
    par.filename = [CSCfiles(it).name];
    if dim(1) > dim(2)% the number of spike events must larger than the number of samples in the waveform
        [sortedSpikes,wltCoeffs] = doSpikeSorting_waveclus( waveclus , par );
    end
    
    % force cluster-membership
    selIx1 = find(sortedSpikes.assignedCluster ==0);
    selIx2 = find(sortedSpikes.assignedCluster ~=0);
    [class_out] = force_membership_wc(wltCoeffs(selIx2,:), sortedSpikes.assignedCluster(selIx2), wltCoeffs(selIx1,:), par);
    sortedSpikes.assignedCluster(selIx1(class_out~=0)) = class_out(class_out~=0);
    
    % save data in fieldtrip-like dummy structure. Note that we can use dum1 to extract the LFP data.
    dum1{it} = [];
    dum1{it}.trial{1} = dataSamplesCSC;
    dum1{it}.time{1} = CSCTime;
    dum1{it}.label = {chanLab};
    dum1{it}.cfg.hdr = hdrCSC;
    dum1{it}.fsample = Fs;
    
    % save the sorted spike-data into fieldtrip-like spike structure
    dum2{it} = [];
    dum2{it}.hdr = hdrCSC;
    dum2{it}.waveform{1}(1,:,:) = sortedSpikes.wavf'; % instead of "spikeWaveforms'"
    dum2{it}.timestamp{1} = SPKTime(sortedSpikes.newSpikeTimes); %% this is weird
    dum2{it}.unit{1} = sortedSpikes.assignedCluster;
    dum2{it}.waveformdimord = '{chan}_lead_time_spike';
    dum2{it}.waveformtime = wft;
    dum2{it}.label = {['sig_',chanLab,'_sorted_wvf']};
    dum2{it}.timestampSec{1} = sortedSpikes.newSpikeTimes/32000; %% I think i added this
    dum2{it}.timestampSample{1} = sortedSpikes.newSpikeTimes;
    
    for ia=1:max(dum2{1,it}.unit{1,1}) % for each single unit in the data set
        idx=find(dum2{1,it}.unit{1,1}==ia);
        spikeTimesCluster=dum2{1,it}.timestampSec{1,1}(idx);
        clusterCount=clusterCount+1;
        % ISI
        expLength=ceil(max(CSCTime));
        dt=linspace(0,expLength,expLength*1000);
        spikeHistogram=hist(spikeTimesCluster,dt);
        if max(spikeHistogram)>1
            disp('Problem with Histogram: Multiple spikes per bin!');
            break
        end
        xc=xcorr(spikeHistogram,250); % maxLag=250ms
        
        close all
        myfig=figure(1);
        numberCluster=max(dum2{1,it}.unit{1,1});
        numberSpikes=sum(dum2{1,it}.unit{1,1}==ia);
        subplot(13,3,[22:3:37],'align');
        hold on
        box on
        bar(xc(252:end), 'b');
        xlabel('Lag [ms]')
        ylabel('Coincidences');
        axis tight
        shortISI=double(sum(xc(252:253))/numberSpikes*100);
        title(sprintf('ISI | Below 3ms: %0.2f%%', shortISI));
        
        
        %% waveform
        waveformCluster=[];
        for ix=1:size(idx,2)
            waveformCluster(ix,:)=dum2{1,it}.waveform{1,1}(:,:,idx(ix));
        end
        
        % plotting waveform
        subplot(13,3,[1:3:16],'align');
        hold on
        for ix=1:size(waveformCluster,1)
            plot(waveformCluster(ix,:),'r');
        end
        plot(mean(waveformCluster,1),'k');
        xlabel('Time [samples]');
        ylabel('Amplitude');
        set(gca,'Ylim',[-100 150]);
        axis tight
        box on
        anzahlCluster=max(dum2{1,it}.unit{1,1});
        title(sprintf('Waveform | Spike %d of %d | N: %d', ia, anzahlCluster, numberSpikes));
        
        %% timewindows
        timeWindow=[];
        timeWindow.encCueLocked=[-1 5];
        timeWindow.encRespLocked=[-2 1];
        timeWindow.retCueLocked=[-1 2];
        timeWindow.retRespLocked=[-2 1];
        
        % add spiketimes of the cluster to the table
        % Encoding + CueLocked
        locking=1; % for CueLocked
        encSpiketimes_cueLocked{size(encSpiketimes_cueLocked,1)+1,:}=insertSpiketimes(encTrigger,spikeTimesCluster,locking,timeWindow.encCueLocked);
        
        % Encoding + RespLocked
        locking=3; % for RespLocked
        encSpiketimes_respLocked{size(encSpiketimes_respLocked,1)+1,:}=insertSpiketimes(encTrigger,spikeTimesCluster,locking,timeWindow.encRespLocked);
        
        % Retrieval + CueLocked
        locking=1;
        retSpiketimes_cueLocked{size(retSpiketimes_cueLocked,1)+1,:}=insertSpiketimes(retTrigger,spikeTimesCluster,locking,timeWindow.retCueLocked);
        
        % Retrieval + RespLocked
        locking=3;
        retSpiketimes_respLocked{size(retSpiketimes_respLocked,1)+1,:}=insertSpiketimes(retTrigger,spikeTimesCluster,locking,timeWindow.retRespLocked);
        
        %% Rasterplots + Firing Frequency
        % timewindows
        timeWindow=[];
        timeWindow.encCueLocked=[-1.5 5.5];
        timeWindow.encRespLocked=[-2.5 1.5];
        timeWindow.retCueLocked=[-1.5 2.5];
        timeWindow.retRespLocked=[-2.5 1.5];
        
        % Rasterplot + Firing frequency: Encoding + CueLocked
        locking=1;
        spiketimes=insertSpiketimes(encTrigger,spikeTimesCluster,locking,timeWindow.encCueLocked);
        [n,dt]=makeRasterplot(timeWindow.encCueLocked, allTrials, spiketimes ,[2:3:11]);
        title(['\fontsize{15}' patientCode ' ' chanLab ' Session: ', sessionNum, sprintf(' Cluster# %d', clusterCount), char(10) '\fontsize{10} Rasterplot: Encoding & Cue-locked']);
        makeFfreq(n,dt,timeWindow.encCueLocked,17) % to run for encCueLocked
        
        % Rasterplot + Firing frequency: Encoding + RespLocked
        locking=3; % for RespLocked
        spiketimes=insertSpiketimes(encTrigger,spikeTimesCluster,locking,timeWindow.encRespLocked);
        [n,dt]=makeRasterplot(timeWindow.encRespLocked,allTrials,spiketimes,[23:3:32]);
        title('Rasterplot: Encoding & Resp-locked');
        makeFfreq(n,dt,timeWindow.encRespLocked,38)
        
        % Rasterplot + Firing frequency: Retrieval + CueLocked
        locking=1;
        spiketimes=insertSpiketimes(retTrigger,spikeTimesCluster,locking,timeWindow.retCueLocked);
        [n,dt]=makeRasterplot(timeWindow.retCueLocked,allTrials,spiketimes,[3:3:12]);
        title('Rasterplot: Retrieval & Cue-locked');
        makeFfreq(n,dt,timeWindow.retCueLocked,18)
        
        % Rasterplot + Firing frequencY: Retrieval + RespLocked
        locking=3;
        spiketimes=insertSpiketimes(retTrigger,spikeTimesCluster,locking,timeWindow.retRespLocked);
        [n,dt]=makeRasterplot(timeWindow.retRespLocked,allTrials,spiketimes,[24:3:33]);
        title('Rasterplot: Retrieval & Resp-locked');
        makeFfreq(n,dt,timeWindow.retRespLocked,39)
        
        % Is what you see a single unit, multiunit or noise?
        set(myfig,'Position',get(0,'Screensize')); % makes it fullscreen
                beep
        

                clusterCategory=''; % make the variable if it doesnt exist yet
                clusterCategory(size(clusterCategory,1)+1,1)=input('Do you see a single unit (s), multiunit (m) or just noise (n)? ','s');
                while clusterCategory(size(clusterCategory,1))~='m' && clusterCategory(size(clusterCategory,1))~='n' && clusterCategory(size(clusterCategory,1))~='s'
                    if size(clusterCategory,1)==1
                        clusterCategory='';  % clear response if was not SMN (first response)
                    else
                        clusterCategory(size(clusterCategory,1))=''; % only clear last response if it was not SMN
                    end
                    disp('Use only the allowed keys (''s'', ''m'', ''n'')!!')
                    clusterCategory(size(clusterCategory,1)+1,1)=input('Do you see a single unit (s), multiunit (m) or just noise (n)? ','s');
                end
%                 save figure and close it
                saveas(myfig, ['allvis_', chanLab, sprintf('_cluster%d',ia),clusterCategory(end), '.png']);
                close all
                
    end
end

save('encSpiketimes_cueLocked','encSpiketimes_cueLocked');
save('encSpiketimes_respLocked','encSpiketimes_respLocked');
save('retSpiketimes_cueLocked','retSpiketimes_cueLocked');
save('retSpiketimes_respLocked','retSpiketimes_respLocked');

%% add SMN to Spiketimes table
encSpiketimes_cueLocked{1,allTrials+1}={'NaN'};
encSpiketimes_cueLocked{2,allTrials+1}={'NaN'};
encSpiketimes_cueLocked.Properties.VariableNames(allTrials+1)={'SMN'};
encSpiketimes_cueLocked{3:end,allTrials+1}=cellstr(clusterCategory);

encSpiketimes_respLocked{1,allTrials+1}={'NaN'};
encSpiketimes_respLocked{2,allTrials+1}={'NaN'};
encSpiketimes_respLocked.Properties.VariableNames(allTrials+1)={'SMN'};
encSpiketimes_respLocked{3:end,allTrials+1}=cellstr(clusterCategory);

retSpiketimes_cueLocked{1,allTrials+1}={'NaN'};
retSpiketimes_cueLocked{2,allTrials+1}={'NaN'};
retSpiketimes_cueLocked.Properties.VariableNames(allTrials+1)={'SMN'};
retSpiketimes_cueLocked{3:end,allTrials+1}=cellstr(clusterCategory);

retSpiketimes_respLocked{1,allTrials+1}={'NaN'};
retSpiketimes_respLocked{2,allTrials+1}={'NaN'};
retSpiketimes_respLocked.Properties.VariableNames(allTrials+1)={'SMN'};
retSpiketimes_respLocked{3:end,allTrials+1}=cellstr(clusterCategory);

% transform spiketimes into number of spikes
% consider preallocating for speed
for ia=1:allTrials
    for ib=3:clusterCount+2
        encSpikeNumber_cueLocked{ib,ia}=num2cell(cellfun(@length,encSpiketimes_cueLocked{ib,ia}));
        encSpikeNumber_respLocked{ib,ia}=num2cell(cellfun(@length,encSpiketimes_respLocked{ib,ia}));
        retSpikeNumber_cueLocked{ib,ia}=num2cell(cellfun(@length,retSpiketimes_cueLocked{ib,ia}));
        retSpikeNumber_respLocked{ib,ia}=num2cell(cellfun(@length,retSpiketimes_respLocked{ib,ia}));
    end
end

encSpikeNumber_cueLocked=encSpikeNumber_cueLocked./2;

% Add SMN to number-table
encSpikeNumber_cueLocked{:,allTrials+1}=encSpiketimes_cueLocked{:,allTrials+1};
encSpikeNumber_cueLocked.Properties.VariableNames(allTrials+1)={'SMN'};

encSpikeNumber_respLocked{:,allTrials+1}=encSpiketimes_respLocked{:,allTrials+1};
encSpikeNumber_respLocked.Properties.VariableNames(allTrials+1)={'SMN'};

retSpikeNumber_cueLocked{:,allTrials+1}=retSpiketimes_cueLocked{:,allTrials+1};
retSpikeNumber_cueLocked.Properties.VariableNames(allTrials+1)={'SMN'};

retSpikeNumber_respLocked{:,allTrials+1}=retSpiketimes_respLocked{:,allTrials+1};
retSpikeNumber_respLocked.Properties.VariableNames(allTrials+1)={'SMN'};

% flag up cluster that need to be ignored
% here multiunits ('m') and noise ('n') is ignored in the next mk_tempCell
% step
myFlag=[];
for i=1:clusterCount
    if strcmp(retSpikeNumber_cueLocked{i+2,allTrials+1},'m') || strcmp(retSpikeNumber_cueLocked{i+2,allTrials+1},'n')
        myFlag(i,1)=1;
    else
        myFlag(i,1)=0;
    end
end

% cueLocked Hits
[tempCell_enc, tempCell_ret]=mk_tempCell(allTrials, encSpikeNumber_cueLocked, retSpikeNumber_cueLocked, 'hit', myFlag);

clusterMean=[];
clusterStd=[];
for i=1:size(tempCell_enc,1)
    clusterMean(i,1)=mean([tempCell_enc{i,:} tempCell_ret{i,:}]);
    clusterStd(i,1)=std([tempCell_enc{i,:} tempCell_ret{i,:}]);
end

RSAhits_cueLocked=normSpikeNumber(tempCell_enc, tempCell_ret, clusterMean, clusterStd);

% cueLocked Misses
[tempCell_enc, tempCell_ret]=mk_tempCell(allTrials, encSpikeNumber_cueLocked, retSpikeNumber_cueLocked, 'miss', myFlag);
RSAmiss_cueLocked=normSpikeNumber(tempCell_enc, tempCell_ret, clusterMean, clusterStd);

% respLocked Hits
[tempCell_enc, tempCell_ret]=mk_tempCell(allTrials, encSpikeNumber_respLocked, retSpikeNumber_respLocked, 'hit', myFlag);
RSAhits_respLocked=normSpikeNumber(tempCell_enc, tempCell_ret, clusterMean, clusterStd);

% respLocked Misses
[tempCell_enc, tempCell_ret]=mk_tempCell(allTrials, encSpikeNumber_respLocked, retSpikeNumber_respLocked, 'miss', myFlag);
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
glmX=zeros(size(hitsIdx,1),3);
for i=1:size(hitsIdx,1)*size(hitsIdx,1)
    glmY(i)=RSAhits_respLocked(i);
    if RSA_mask(i)==1
        glmX(i,1)=1;
    elseif RSA_mask(i)==2
        glmX(i,2)=1;
    elseif RSA_mask(i)==3
        glmX(i,3)=1;
    end
end

RSAcoeff_hitsResp=glmfit(glmX,glmY,'normal');

% GLM for hits + cueLocked
glmX=zeros(size(hitsIdx,1),3);
glmY=[];
for i=1:size(hitsIdx,1)*size(hitsIdx,1)
    glmY(i)=RSAhits_cueLocked(i);
    if RSA_mask(i)==1
        glmX(i,1)=1;
    elseif RSA_mask(i)==2
        glmX(i,2)=1;
    elseif RSA_mask(i)==3
        glmX(i,3)=1;
    end
end

RSAcoeff_hitsCue=glmfit(glmX,glmY,'normal');

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

glmX=zeros(size(missIdx,1),3);
glmY=[];
for i=1:size(missIdx,1)*size(missIdx,1)
    glmY(i)=RSAmiss_cueLocked(i);
    if RSA_mask(i)==1
        glmX(i,1)=1;
    elseif RSA_mask(i)==2
        glmX(i,2)=1;
    elseif RSA_mask(i)==3
        glmX(i,3)=1;
    end
end

RSAcoeff_missCue=glmfit(glmX,glmY,'normal');

glmX=zeros(size(missIdx,1),3);
glmY=[];
for i=1:size(missIdx,1)*size(missIdx,1)
    glmY(i)=RSAmiss_respLocked(i);
    if RSA_mask(i)==1
        glmX(i,1)=1;
    elseif RSA_mask(i)==2
        glmX(i,2)=1;
    elseif RSA_mask(i)==3
        glmX(i,3)=1;
    end
end

RSAcoeff_missResp=glmfit(glmX,glmY,'normal');