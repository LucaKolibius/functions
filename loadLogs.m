function [tableTemplate, hitsIdx, missIdx, allTrials, sessionNum, retTrigger, encTrigger, orda, animalCues] = loadLogs(p2d, trls)
% error('ERLANGEN TRIGGER ARE SAMPLED AT ~32678 Hz');
% load in TTLs
cd(p2d)
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

% [TimeStampsTTL, ttls, HdrTTL] = Nlx2MatEV(EVfile.name, FieldSelection, ExtractHeader, ExtractMode, ModeArray );
cd('\\analyse4.psy.gla.ac.uk\project0309\Luca\toolboxes\Neuralynx_19012019\') % after restarting the neuralynx function does not find the path to a mex file, unless you cd into that folder
[TimeStampsTTL, ttls, ~] = Nlx2MatEV([EVfile.folder, filesep, EVfile.name], FieldSelection, ExtractHeader, ExtractMode, ModeArray );
% [timestampsCSC, dataSamplesCSC,hdrCSC] = Nlx2MatCSC([p2d,CSCfiles(it).name], FieldSelection, ExtractHeader, ExtractMode, []);% import the raw data

cd(p2d)

events = zeros(size(ttls,2),2);
events(:,1) = TimeStampsTTL';
events(:,2) = ttls';
TimeStampsTTL = TimeStampsTTL-TimeStampsTTL(1,1);
%TimeStampsTTL = TimeStampsTTL./31.25; % convert from us to samples
TimeStampsTTL = TimeStampsTTL./1e6; % convert from us to seconds
events(:,1) = TimeStampsTTL'; % update event variable with new times (which are now in seconds starting at 0)

% % deletes 7s that appear within a 255 or 128 TTL
% last128 = find(events(:,2)==128);
% last255 = find(events(:,2)==255);
% 
% % i don't know why the if statement is suddenly necessary when there are no
% % 128's or 255's
% if ~isempty(last128) || ~isempty(last255)
%     lastInit = max([last128(end), last255(end)]); % finds the last index at which a 128 or 225 TTL is sent
% else
%     lastInit= [];
% end

a = find(events(:,2)==7);
% i used this when we had another session in the same
% recording and we had 128 and 255 TTl pulses afterwards because we tested
% the trigger after a crashed experiment session anew
% a(62:end)=[]; 

% % exclude this for P08ERL. Somehow 128 and 255 trigger are at the end of
% % the recording.
% if isempty(regexp(p2d,'P08ERL','ONCE')) % if it is not P08ERL, continue with L38, otherwise skip this part
%     if isempty(lastInit)
%         warning('No 128 or 255 TTL pulses');
%     else
%         a(a<=lastInit)=[]; % if there are no 128 or 255 TTL pulses, this line will result in an error
%     end
% end

% % extract timestamp of the trial start
trialStart = events(a,1); % on which index in the variable "events" do we have the trigger #7

% only for P05_S4 (I had two sessions within one recording)
% this setting (trialStart = trialStarta) is for S4; S5 does not yet
% reconstruct properly (15.7.19)
% if ~isempty(regexp(p2d,'P05', 'ONCE')) && ~isempty(regexp(p2d,'S4', 'ONCE')) && isempty(regexp(p2d,'ERL', 'ONCE')) % P05_S4, not P05ERL_S4
%     % old from when two sessions were in the same recording
%     trialStarta = trialStart(1:98);
%     trialStartb = trialStart(99:196);
%     trialStart=trialStarta;
% end

%% logfile
% dir('*fVSp*txt'); % shows all logfiles in current directory
% logFilename = input('What is the name of the Log-File? ', 's');

mDir = dir('*fVSp*txt');
logFilename = mDir.name;

patientCode = logFilename(1:6); %for ERL suffix
patientCode(regexp(patientCode,'_fV'):regexp(patientCode,'_fV')+2) = []; % delete part if there was no ERL suffix
startSessionNum = regexp(logFilename,'fVSp_')+5;
endSessionNum = regexp(logFilename,'_');
endSessionNum = endSessionNum(4)-1;
sessionNum = logFilename(startSessionNum:endSessionNum); % not all logfiles have the same grammar, so this should prevent problems
allLogfiles = dir('*fVSp*txt');
numLogfiles = size(allLogfiles, 1);
disp([patientCode, sessionNum]);

if numLogfiles ~= 1
    error('not sure if the xcorr works for a spurious ttl or logfile entry withou ttl')
end

sessionNum(regexp(sessionNum,'_'))='-'; % changes '_' to a '-' so it does not mess up my title

% read in logfile
raw = [];
for ixx=numLogfiles:-1:1 % this is not 1:numLogfiles because for whatever reason in P07ERL_S2 the later logfile is above the first logfile. I circumvent this by just reversing the call order in the following loop.
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

loglist = cellfun(@str2double, raw_temp(:,1)); % new
% % finds out which trials are presented during encoding AND retrieval, saves this
% % information as a binary code into log_trialSave
log_trialSave = zeros(size(loglist,1),1);
for i=1:size(loglist,1) % numLog/2 = trials according to logfile
    if size(find(loglist == loglist(i)),1) == 2 % if the trial # appears two times in the logfile
        log_trialSave(i) = 1;
    elseif isnan(loglist(i))
        log_trialSave(i) = 1;
    end
end

% if ~isempty(regexp(p2d,'P09ERL','ONCE')) || ( ~isempty(regexp(p2d,'P07ERL','ONCE')) && ~isempty(regexp(p2d,'S1','ONCE')) ) % for P09ERL and P07ERL_S1
%     datalist=loglist;
%     datalist(isnan(datalist))=[];
%     data_trialSave = zeros(size(trialStart,1),1);
%     for i=1:size(datalist,1) % numLog/2 = trials according to logfile
%         if size(find(datalist==datalist(i)),1) == 2
%             data_trialSave(i)=1;
%         end
%     end
% 
% else
% CHECK IF TRUE ALL THE TIME!!
% new implementation. here the default is to keep the ttl in, unless there
% the trial appears twice (once in encoding and once in retrieval). because
% in an incomplete logfile most trials will not have a retrieval trial
% datalist=loglist;
% datalist(isnan(datalist))=[];
% data_trialSave = ones(size(trialStart,1),1);
% for i=1:size(datalist,1) % numLog/2 = trials according to logfile
%     if size(find(datalist==datalist(i)),1) ~= 2
%         data_trialSave(i)=0;
%     end
% end
% end

% % loglist is with nan's; datalist is basically triggers infered from logfile, without nans
% if ~isempty(regexp(p2d,'P05ERL','ONCE')) && ~isempty(regexp(p2d,'S1','ONCE')) % fix for P05ERLS1 (there was no trigger for the first trial)
%     delLog = find(loglist==1); % index of first trial in enc and ret
%     log_trialSave(delLog) = 0; % delete trial 1 encoding and retrieval from the logfile
%     
%     delTrig = find(datalist==1);
%     delTrig = delTrig(2)-1; % the first trigger got lost. because of that the corresponding trigger for retrieval-trial-1 is at position 30 instead of 31
%     data_trialSave(delTrig) = 0;
% end
    
% delete the entries in the logfiles that cannot be saved according to
% log_trialSave:
disp(sprintf('Deleting %d logfiles', size(raw_temp(log_trialSave==0,:),1) ));
% disp(sprintf('Deleting %d logfiles and %d TTL trigger', size(raw_temp(log_trialSave==0,:),1), size(trialStart(data_trialSave==0),1) ));
raw_temp(log_trialSave==0,:)=[];
% trialStart(data_trialSave==0)=[];
%%%

raw = [raw; raw_temp];
end

% numTTL = size(trialStart,1);

% old
% numLog = max(str2double(raw_temp(:,1)) - length(find(log_trialSave==0)) )*2; % find the maximum number and multiply by two (one for encoding and retrieval each)

% if ~isempty(regexp(p2d,'P05ERL','ONCE')) && ~isempty(regexp(p2d,'S1','ONCE'))
%     numLog = numTTL; % numLog is not properly calculated because in P05ERL_S1 I take out a trial from the beginning. Calculating the number of logfile entries through the max number * 2 does not work because of that. This is a lazy fix.
% else
%     numLog = max(str2double(raw(:,1))*2); % find the maximum number and multiply by two (one for encoding and retrieval each)
% end

% there is more data than logfile entries. this implies that there was a 
% previoius logfile with faulty data. therefore we prune the data files 
% to this extent. This approach does not work anymore because apparently 
% some logfiles still represent valid trials. Ergo I will first extract 
% these valid trials, delete the rest from the logfile and the data file 
% (so trialStart) and add the size of the previous logfile (session 1a) 
% to the next logfile (session 1b), so they are in sync with the data 
% (trialStart)

% % if the number of TTL spikes are not equal to the trigger in the log file,
% if numTTL~=numLog
%     disp('Dude! Unterschiedliche Anzahl von TTL und Logfile Trigger');
%     % find direction of difference
%     if numTTL>numLog
%         disp('There are more TTL spikes than Logfile entries. This shouldn''t be - time to worry! (break at line l42)');
%         %         diffSize=numTTL-numLog;
%         %         trialStart(1:diffSize)=[];
%         return
%     elseif numTTL<numLog
%         disp('There are more Logfile entries than TTL spikes. This is the time to worry! (break at line 147)');
%         return
%     end
% end

%% Separate encoding and retrieval
% delLog=str2double(raw(:,1)); % old; doesn't work anymore for whatever reason
delLog = cellfun(@str2double, raw(:,1)); % new


while ~isnan(delLog(end-1)) % some logfiles (P09S03) do not have NaN values at the end. this is added here because this script relies on them
    delLog(end+1)=NaN;
end

flag=0; % flag 0/1 -> encoding/retrieval
blockcounter=1;
idxEnc=[];
idxRet=[];
for ix=1:size(delLog,1)
    if ix==1
        encBlock_start=1;
    elseif isnan(delLog(ix)) && isnan(delLog(ix-1))
        encBlock_start=ix+1;
    elseif isnan(delLog(ix)) && isnumeric(delLog(ix+1)) && flag==0 % encoding blocks
        idxEnc=[idxEnc; [encBlock_start:ix-1]'];
        flag=1;
        retBlock_start=ix+1;
    elseif isnan(delLog(ix)) && flag==1 % retrieval blocks
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

%% determine trialStartLog
trialStartEnc = cellfun(@str2double, [rawEnc(:,5)]);
trialStartRet = cellfun(@str2double, [rawRet(:,10)]);
trialStartLog = sort([trialStartEnc; trialStartRet]);

%% determine ttl and logfile match
preSizeTTL = length(trialStart);
preSizeLog = length(trialStartLog);
if ~isequal(length(trialStartLog), length(trialStart))
    %% visu - pre
%     figure
%     subplot(211)
%     dt = linspace(trialStartLog(1), trialStartLog(end), 10000);
%     histogram(trialStartLog,dt);
%     subplot(212)
%     dt = linspace(trialStart(1), trialStart(end), 10000);
%     histogram(trialStart,dt);
    
    
    %%
    [trgIdx, stopIt] = trigger_convolution(trialStartLog, trialStart);
    if stopIt == -1 % more ttl than logfiles
    trialStart = trialStart(trgIdx : trgIdx+size(trialStartLog,1)-1);
    elseif stopIt == 1 % more logfiles than ttl
        trialStartLog = trialStartLog(trgIdx : trgIdx+size(trialStart,1)-1); % prune logfiles trigger
        delTag = [raw{1:trgIdx-1}]; % tag the logfile that does not have a TTL, so later we can get rid of the retrieval logfile and trigger
        raw(1:trgIdx-1,:) = []; % cut out the logfile entry that does not have a TTL from "raw"
        
        if trgIdx ~= 1 % if the logfiles at the beginning of trialStartLog need to be pruned, the corresponding retrieval logfiles and TTLs have to be eliminated
            delLogs = ismember([raw{:,1}],delTag)'; % index of where in the raw-file the corresponding retrieval trial sits
            delStarts = delLogs; % the index for the trialStart and trialStartLogs is almost the same ...
            delStarts(isnan(str2double(raw(:,1)))) = []; % ... we just have to get rid of the nans that we have in the raw file
            
            % wipe out the correpsonding retrieval trials from the
            % rawfile/trialStartLog/trialStart
            raw(delLogs,:) = [];  % delete the corresponding retrieval trial from the raw file     
            trialStartLog(delStarts) = [];
            trialStart(delStarts) = [];
        end
%     if stopIt == 1
%         error('More Logfiles than triggers! Investigate')
%     end
    end
     %% visu - post
%     figure
%     subplot(211)
%     dt = linspace(trialStartLog(1), trialStartLog(end), 10000);
%     histogram(trialStartLog,dt);
%     subplot(212)
%     dt = linspace(trialStart(1), trialStart(end), 10000);
%     histogram(trialStart,dt);
end

postSizeTTL = length(trialStart);
postSizeLog = length(trialStartLog);
disp(sprintf('Deleted %d TTLs', preSizeTTL-postSizeTTL));
disp(sprintf('Deleted %d Logfile entries', preSizeLog-postSizeLog));


%% add trialStart (in s) to raw logFile
counter=1; % using a counter here to avoid problems with NaN
for ix=1:size(raw,1)
    if ~isnan(str2double(raw{ix,1}))
        raw{ix,13}=trialStart(counter);
        counter=counter+1;
    end
end

%% stimulus place or face?
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

%% separate encoding and retrieval again, now with the trialStart added
% delLog=str2double(raw(:,1)); % old; doesn't work anymore for whatever reason
delLog = cellfun(@str2double, raw(:,1)); % new

while ~isnan(delLog(end-1)) % some logfiles (P09S03) do not have NaN values at the end. this is added here because this script relies on them
    delLog(end+1)=NaN;
end

flag=0; % flag 0/1 -> encoding/retrieval
blockcounter=1;
idxEnc=[];
idxRet=[];
for ix=1:size(delLog,1)
    if ix==1
        encBlock_start=1;
    elseif isnan(delLog(ix)) && isnan(delLog(ix-1))
        encBlock_start=ix+1;
    elseif isnan(delLog(ix)) && isnumeric(delLog(ix+1)) && flag==0 % encoding blocks
        idxEnc=[idxEnc; [encBlock_start:ix-1]'];
        flag=1;
        retBlock_start=ix+1;
    elseif isnan(delLog(ix)) && flag==1 % retrieval blocks
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

%% Trigger for encoding and retrieval (1: cue; 2: stimulus; 3: response)
encTrigger = cellfun(@str2double, rawEnc(:,5:7)); % new

for ix=1:size(rawEnc,1)
    subst=encTrigger(ix,1)-cell2mat(rawEnc(ix,13));
    encTrigger(ix,1)=encTrigger(ix,1)-subst;
    encTrigger(ix,2)=encTrigger(ix,2)-subst;
    encTrigger(ix,3)=encTrigger(ix,3)-subst;
end

retTrigger = cellfun(@str2double, rawRet(:,10:12)); % new
for ix=1:size(rawEnc,1)
    subst=retTrigger(ix,1)-cell2mat(rawRet(ix,13));
    retTrigger(ix,1)=retTrigger(ix,1)-subst;
    retTrigger(ix,2)=retTrigger(ix,2)-subst;
    retTrigger(ix,3)=retTrigger(ix,3)-subst;
end
allTrials=size(encTrigger,1);

% enc1 und ret1 should be the same stimuli
retTrigger=[retTrigger cellfun(@str2double, rawRet(:,1))];
retTrigger=sortrows(retTrigger,4);
retTrigger = retTrigger(1:end,1:3); % delete the sorting column

% hits or misses
% find out trials with two hits
respCorr = cellfun(@str2double, rawRet(:,[1 5 6]));
respCorr = sortrows(respCorr,1);

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

%% making the table
tempCell=[];
trials=[];
for i=1:allTrials
    trials{1,i}=sprintf('trial%d',i);
    value=[faceORplace(i,1),faceORplace(i,2)];
    value=strjoin(value,''); % combine the two characters, use no delimiter
    tempCell{1,i}=value;
    tempCell{2,i}=faceORplace(i,3);
end
tableTemplate = cell2table(tempCell);
tableTemplate.Properties.VariableNames = trials(1,:);

% this is an index to reorder retrieval according to its presentation order
orda = cellfun(@str2double, rawRet(:,1));

idx = ~cellfun(@isempty, strfind(tableTemplate{2,:}, 'miss'));
orda(idx,:) = [];

orda(:,2) = [1:size(orda,1)]';
orda = sortrows(orda,1);
orda = orda(:,2);

%% animal cue
% trls is refering to hits only
animalCues = [rawEnc{hitsIdx(trls),2}];
end

