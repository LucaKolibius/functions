function [mCell, allTrials, trlTrigger] = loadLogs_tEMt(p2d)
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

[TimeStampsTTL, ttls, ~] = Nlx2MatEV(EVfile.name, FieldSelection, ExtractHeader, ExtractMode, ModeArray );
events = zeros(size(ttls,2),2);
events(:,1) = TimeStampsTTL';
events(:,2) = ttls';
TimeStampsTTL = TimeStampsTTL-TimeStampsTTL(1,1);
%TimeStampsTTL = TimeStampsTTL./31.25; % convert from us to samples
TimeStampsTTL = TimeStampsTTL./1e6; % convert from us to seconds
events(:,1) = TimeStampsTTL'; % update event variable with new times (which are now in seconds starting at 0)

% deletes 7s that appear within a 255 or 128 TTL
last128 = find(events(:,2)==128);
last255 = find(events(:,2)==255);

% i don't know why the if statement is suddenly necessary when there are no
% 128's or 255's
if ~isempty(last128) || ~isempty(last255)
    lastInit = max([last128(end), last255(end)]); % finds the last index at which a 128 or 225 TTL is sent
else
    lastInit= [];
end

a = find(events(:,2)==10);
% i used this when we had another session in the same
% recording and we had 128 and 255 TTl pulses afterwards because we tested
% the trigger after a crashed experiment session anew
% a(62:end)=[]; 

if isempty(lastInit)
    warning('No 128 or 255 TTL pulses');
else
    a(a<=lastInit)=[]; % if there are no 128 or 255 TTL pulses, this line will result in an error
end

% % extract timestamp of the trial start
trialStart = events(a,1); % on which index in the variable "events" do we have the trigger #7

%% logfile
% dir('*fVSp*txt'); % shows all logfiles in current directory
% logFilename = input('What is the name of the Log-File? ', 's');

mDir = dir('*ctune*txt');
logFilename = mDir.name;
allLogfiles = dir('*ctune*txt');
numLogfiles = size(allLogfiles, 1);

raw = [];
for ixx=numLogfiles:-1:1 % this is not 1:numLogfiles because for whatever reason in P07ERL_S2 the later logfile is above the first logfile. I circumvent this by just reversing the call order in the following loop.

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

% %% NEW STUFF
% % Sometimes incomplete log files have usable trials. This new section
% % should harvest that data
% 
% loglist = str2double(raw_temp(:,1));
% % finds out which trials are presented during encoding AND retrieval, saves this
% % information as a binary code into log_trialSave
% log_trialSave = zeros(size(loglist,1),1);
% for i=1:size(loglist,1) % numLog/2 = trials according to logfile
%     if size(find(loglist == loglist(i)),1) == 2 % if the trial # appears two times in the logfile
%         log_trialSave(i) = 1;
%     elseif isnan(loglist(i))
%         log_trialSave(i) = 1;
%     end
% end
% 
% % new implementation. here the default is to keep the ttl in, unless there
% % the trial appears twice (once in encoding and once in retrieval). because
% % in an incomplete logfile most trials will not have a retrieval trial
% datalist=loglist;
% datalist(isnan(datalist))=[];
% data_trialSave = ones(size(trialStart,1),1);
% for i=1:size(datalist,1) % numLog/2 = trials according to logfile
%     if size(find(datalist==datalist(i)),1) ~= 2
%         data_trialSave(i)=0;
%     end
% end
% 
% % delete the entries in the logfiles that cannot be saved according to
% % log_trialSave:
% disp(sprintf('Deleting %d logfiles and %d TTL trigger', size(raw_temp(log_trialSave==0,:),1), size(trialStart(data_trialSave==0),1) ));
% raw_temp(log_trialSave==0,:)=[];
% trialStart(data_trialSave==0)=[];
% %%%

raw = [raw; raw_temp];
end

numTTL = size(trialStart,1);
numLog = size(raw,1); % find the maximum number and multiply by two (one for encoding and retrieval each)

% if the number of TTL spikes are not equal to the trigger in the log file,
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

% add trialStart (in s) to raw logFile
counter=1; % using a counter here to avoid problems with NaN
for ix=1:size(raw,1)
    if ~isnan(str2double(raw{ix,1}))
        raw{ix,13}=trialStart(counter);
        counter=counter+1;
    end
end

allTrials = size(raw,1);
trlTrigger = zeros(allTrials,1);
for ib = 1:size(raw,1)
    subst = str2double(raw(ib,5)) - cell2mat(raw(ib,13));
    trlTrigger(ib,1) = str2double(raw(ib,6)) - subst; % first flip of image
end

% %% making the table
% tempCell=[];
% trials=[];
% for i=1:allTrials
%     trials{1,i}=sprintf('trial%d',i);
%     value=[faceORplace(i,1),faceORplace(i,2)];
%     value=strjoin(value,''); % combine the two characters, use no delimiter
%     tempCell{1,i}=value;
%     tempCell{2,i}=faceORplace(i,3);
% end
% tableTemplate = cell2table(tempCell);
% tableTemplate.Properties.VariableNames = trials(1,:);

%% 
raw(:,14) = {0};
imgName = {};
timestamps = {};
mCell = {};
for ix = 1:size(raw,1)
  
    if raw{ix,14} == 1
        continue
    end
    
    mCell(1,end+1) = cellstr(raw{ix,2});
    mCell(2,end)   = {find(strcmp([raw{:,2}],mCell(1,end)))};
    raw(cell2mat(mCell(2,end)), 14) = {1};
end

end % end of function