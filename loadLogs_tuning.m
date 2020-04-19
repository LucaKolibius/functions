function [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, hitsIdx, missIdx, allTrials, sessionNum, retTrigger, encTrigger]=loadLogs(p2d)
% retTrigger and encTrigger as output?!

%trigger
FieldSelection(1) = 1; %timestamps
FieldSelection(2) = 0; % EventIDs
FieldSelection(3) = 1; %TTLs
FieldSelection(4) = 0; % Extras
FieldSelection(5) = 0; % Event strings
ExtractHeader = 1;
ExtractMode = 1;
ModeArray = [];
EVfile = dir([p2d,'*nev']);

[TimeStampsTTL, ttls, HdrTTL] = Nlx2MatEV(EVfile.name, FieldSelection, ExtractHeader, ExtractMode, ModeArray );
events = zeros(size(ttls,2),2);
events(:,1) = TimeStampsTTL';
events(:,2) = ttls';
TimeStampsTTL = TimeStampsTTL-TimeStampsTTL(1,1);
%TimeStampsTTL = TimeStampsTTL./31.25; % convert from us to samples
TimeStampsTTL = TimeStampsTTL./1e6; % convert from us to seconds
events(:,1) = TimeStampsTTL'; % update event variable with new times (which are now in seconds starting at 0)
% TTLidx = find(events(:,2) == 7, 1); % on which index in the variable "events" do we have the trigger #7
% 
% % extract timestamp of the trial start
% for i=1:length(TTLidx)
%     trialStart(i)=events(TTLidx(i),1);
% end
trialStart = events(events(:,2)==7,1);

%% logfile
dir('*fVSp*txt'); % shows all logfiles in current directory
logFilename = input('What is the name of the Log-File? ', 's');

patientCode = logFilename(1:6); %for ERL suffix
patientCode(regexp(patientCode,'_fV'):regexp(patientCode,'_fV')+2) = []; % delete part if there was no ERL suffix
startSessionNum = regexp(logFilename,'fVSp_')+5;
endSessionNum = regexp(logFilename,'_');
endSessionNum = endSessionNum(4)-1;
sessionNum = logFilename(startSessionNum:endSessionNum); % not all logfiles have the same grammar, so this should prevent problems

allLogfiles = dir('*fVSp*txt');
numLogfiles = size(allLogfiles, 1);

sessionNum(regexp(sessionNum,'_'))='-'; % changes '_' to a '-' so it does not mess up my title

raw = [];
for ixx=1:numLogfiles
% this script imports the logfile into a nice table
% Initialize variables.

filename = allLogfiles(ixx).name;
delimiter = '\t';
startRow = 7;

% Read columns of data as text:
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

% Open the text file.
fileID = fopen(filename, 'r');

% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

% Close the text file.
fclose(fileID);

% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw_temp = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw_temp(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end

%% NEW STUFF
% Sometimes incomplete log files have usable trials. This new section
% should harvest that data

loglist=str2double(raw_temp(:,1));
% finds out which trials are presented during encoding AND retrieval, saves this
% information as a binary code into log_trialSave
log_trialSave=zeros(size(loglist,1),1);
for i=1:size(loglist,1) % numLog/2 = trials according to logfile
    if size(find(loglist==loglist(i)),1)==2
        log_trialSave(i)=1;
    elseif isnan(loglist(i))
        log_trialSave(i)=1;
    end
end

datalist=loglist;
datalist(isnan(datalist))=[];
data_trialSave=zeros(size(trialStart,1),1);
for i=1:size(datalist,1) % numLog/2 = trials according to logfile
    if size(find(datalist==datalist(i)),1)==2
        data_trialSave(i)=1;
    end
end

% delete the entries in the logfiles that cannot be saved according to
% log_trialSave:
disp(sprintf('Deleting %d logfiles and %d TTL trigger', size(raw_temp(log_trialSave==0,:),1), size(trialStart(data_trialSave==0),1)));
raw_temp(log_trialSave==0,:)=[];
trialStart(data_trialSave==0)=[];
%%%

raw=[raw; raw_temp];
end

numTTL=size(trialStart,1);
numLog=max(str2double(raw_temp(:,1)) - length(find(log_trialSave==0)))*2; % find the maximum number and multiply by two (one for encoding and retrieval each)

% there is more data than logfile entries. this implies that there was a 
% previoius logfile with faulty data. therefore we prune the data files 
% to this extent. This approach does not work anymore because apparently 
% some logfiles still represent valid trials. Ergo I will first extract 
% these valid trials, delete the rest from the logfile and the data file 
% (so trialStart) and add the size of the previous logfile (session 1a) 
% to the next logfile (session 1b), so they are in sync with the data 
% (trialStart)

% if the number of TTL spikes are not equal to the trigger in the log file,
if numTTL~=numLog
    disp('Dude! Unterschiedliche Anzahl von TTL und Logfile Trigger');
    % find direction of difference
    if numTTL>numLog
        disp('There are more TTL spikes than Logfile entries. This shouldn''t be - time to worry! (break at line l14)');
        %         diffSize=numTTL-numLog;
        %         trialStart(1:diffSize)=[];
        return
    elseif numTTL<numLog
        disp('There are more Logfile entries than TTL spikes. This is the time to worry! (break at line 129)');
        return
    end
end

%% add trialStart (in s) to raw logFile
counter=1; % using a counter here to avoid problems with NaN
for ix=1:size(raw,1)
    if ~isnan(str2double(raw{ix,1}))
        raw{ix,13}=trialStart(counter);
        counter=counter+1;
    end
end

%% Separate encoding and retrieval
temp=str2double(raw(:,1));
flag=0; % flag 0/1 -> encoding/retrieval
blockcounter=1;
idxEnc=[];
idxRet=[];
for ix=1:size(temp,1)
    if ix==1
        encBlock_start=1;
    elseif isnan(temp(ix)) && isnan(temp(ix-1))
        encBlock_start=ix+1;
    elseif isnan(temp(ix)) && isnumeric(temp(ix+1)) && flag==0 % encoding blocks
        idxEnc=[idxEnc; [encBlock_start:ix-1]'];
        flag=1;
        retBlock_start=ix+1;
    elseif isnan(temp(ix)) && flag==1 % retrieval blocks
        idxRet=[idxRet; [retBlock_start:ix-1]'];
        flag=0;
        blockcounter=blockcounter+1;
    end
end

rawEnc={};
rawRet={};
for ix=1:size(raw,1)
    if find(idxEnc==ix)
        rawEnc=[rawEnc; raw(ix,:)];
    elseif find(idxRet==ix)
        rawRet=[rawRet; raw(ix,:)];
    end
end

% stimulus place or face?
temp1=cellstr(rawEnc(:,3));
stim1={};
for ix=1:size(temp1,1)
    stim1{ix}=temp1{ix,1}(1);
end
stim1=stim1';

temp2=cellstr(rawEnc(:,4));
stim2={};
for ix=1:size(temp2,1)
    stim2{ix}=temp2{ix,1}(1);
end
stim2=stim2';
faceORplace=[stim1 stim2];


%% Trigger for encoding and retrieval (1: cue; 2: stimulus; 3: response)
encTrigger=str2double(rawEnc(:,5:7));
for ix=1:size(rawEnc,1)
    subst=encTrigger(ix,1)-cell2mat(rawEnc(ix,13));
    encTrigger(ix,1)=encTrigger(ix,1)-subst;
    encTrigger(ix,2)=encTrigger(ix,2)-subst;
    encTrigger(ix,3)=encTrigger(ix,3)-subst;
end

retTrigger=str2double(rawRet(:,10:12));
for ix=1:size(rawEnc,1)
    subst=retTrigger(ix,1)-cell2mat(rawRet(ix,13));
    retTrigger(ix,1)=retTrigger(ix,1)-subst;
    retTrigger(ix,2)=retTrigger(ix,2)-subst;
    retTrigger(ix,3)=retTrigger(ix,3)-subst;
end
allTrials=size(encTrigger,1);

% enc1 und ret1 should be the same stimuli
retTrigger=[retTrigger str2double(rawRet(:,1))];
retTrigger=sortrows(retTrigger,4);
retTrigger = retTrigger(1:end,1:3); % delete the sorting column

% hits or misses
% find out trials with two hits
respCorr=str2double(rawRet(:,[1 5 6]));
respCorr=sortrows(respCorr,1);

hitsIdx=[];
for ia=1:size(respCorr,1)
    if respCorr(ia,2) == 1 && respCorr(ia,3) == 1
        hitsIdx(size(hitsIdx,1)+1,1)=ia;
        faceORplace{ia,3}='hit';
    end
end

missIdx=[];
for ia=1:size(respCorr,1)
    if respCorr(ia,2) ~=1 || respCorr(ia,3) ~=1
        missIdx(size(missIdx,1)+1,1)=ia;
        faceORplace{ia,3}='miss';
    end
end

% % trialChange - hits and misses
% % face and face
% temp=[];
% counter=0;
% for ia=1:length(faceORplace)
%     if strcmp(faceORplace(ia,1),'f') && strcmp(faceORplace(ia,2),'f')
%         counter=counter+1;
%         temp(ia,1)=counter;
%     end
% end
% 
% % place and place
% for ia=1:length(faceORplace)
%     if strcmp(faceORplace(ia,1),'p') && strcmp(faceORplace(ia,2),'p')
%         counter=counter+1;
%         temp(ia,1)=counter;
%     end
% end
% 
% % face and place
% for ia=1:length(faceORplace)
%     if strcmp(faceORplace(ia,1),'f') && strcmp(faceORplace(ia,2),'p')
%         counter=counter+1;
%         temp(ia,1)=counter;
%     end
% end
% 
% for ia=1:length(temp)
%     trialChange(ia,1)=find(temp==ia);
% end
% 
% % hits after changetrial
% hitsIdx_chTrial=trialChange(hitsIdx); % these are the hit trials after the changeTrial order switch
% missIdx_chTrial=trialChange(missIdx);

%% making the table
% Encoding + CueLocked
tempCell=[];
trials=[];
for i=1:allTrials
    trials{1,i}=sprintf('trial%d',i);
    value=[faceORplace(i,1),faceORplace(i,2)];
    value=strjoin(value,''); % combine the two characters, use no delimiter
    tempCell{1,i}=value;
    tempCell{2,i}=faceORplace(i,3);
end
encSpiketimes_cueLocked=cell2table(tempCell);
encSpiketimes_cueLocked.Properties.VariableNames=trials(1,:);

% Encoding + RespLocked
encSpiketimes_respLocked=cell2table(tempCell);
encSpiketimes_respLocked.Properties.VariableNames=trials(1,:);

% Retrieval + CueLocked
retSpiketimes_cueLocked=cell2table(tempCell);
retSpiketimes_cueLocked.Properties.VariableNames=trials(1,:);

% Retrieval + RespLocked
retSpiketimes_respLocked=cell2table(tempCell);
retSpiketimes_respLocked.Properties.VariableNames=trials(1,:);

% % pretrial baseline
% encSpikeNumbers_pretrial=double.empty;
% retSpikeNumbers_pretrial=double.empty;

end