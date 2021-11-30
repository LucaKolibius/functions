% find out trial length

encLength = [];
retLength = [];
for sbj = 1:size(allSubj,1)

try
    cd X:/Luca/data
catch
    cd /media/ldk898/rds-share/Luca/data
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
abc = dir;
cd(abc(3).name)
p2d = cd;
p2d(end+1)='\';
[~, hitsIdx, missIdx, allTrials, sessionNum, retTrigger, encTrigger] = loadLogs(p2d);

subjID = allSubj{sbj};
for ia = 1:size(retTrigger,1)
    encLength(1, size(encLength,2)+1) = encTrigger(ia,3) - encTrigger(ia,1);
    retLength(1, size(retLength,2)+1) = retTrigger(ia,3) - retTrigger(ia,1);
end
end

x = 1:1:120;
figure(1)
hist(encLength,x)
min (encLength)
max (encLength)
ylim([0 350]);

figure(2)
hist(retLength,x)
min (retLength)
max (retLength)

%% 40s limit % 500ms TW
% manualRej -> clNames
% postXCrej -> tableTimestamps_posstXC
try
    cd X:/Luca/data
catch
    cd /media/ldk898/rds-share/Luca/data
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
abc = dir;
cd(abc(3).name)
cd advancedAnalysis\manualRej\
abc = dir('clNames*');
load(abc.name)

cd ../postXCrej
abc = dir('tableTimestamps_postXC*');
load(abc.name)

twidth = 0.5; % 500ms
endTW = 41/twidth; % 41s (-1 to cue until 40s after)
for trl = 1:size(encTrigger,1) % trials
    if encTrigger(trl,3)-encTrigger(trl,1) > 40 || retTrigger(trl,3)-retTrigger(trl,1) > 40 % if encoding or retrieval trials take more than 40 seconds, continue
        continue
    end   
    
    % if the time window goes into next trial, continue
    if any(extr)>encTrigger(
    
    for tw = 1:endTW
        timeWindow = -1:twidth:40;
        
        extr = [timeWindow(tw) timeWindow(tw+1)];

        temp = insertSpiketimes(encTrigger, x, 1, extr); % insertspikes distributes the spikes from ONE single units into trials
    end
ENC.TW500ms(1,2).trial = zeros(5,5);