clear all
p2d = cd;
p2d(end+1) = '\';
par = set_parameters;
CSCfiles = dir([p2d,'*.ncs']);

it = 63;

% define the FieldSelection variable and set the ExtracMode.
FieldSelection(1) = 1;%timestamps
FieldSelection(2) = 0;
FieldSelection(3) = 0;%sample freq
FieldSelection(4) = 0;
FieldSelection(5) = 1;%samples
ExtractHeader = 1;
ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.


[timestampsCSC, dataSamplesCSC, hdrCSC] = Nlx2MatCSC(CSCfiles(it).name, FieldSelection, ExtractHeader, ExtractMode, []);
% extract the scale factor
chck = regexp(hdrCSC,'ADBitVolts');
selIdx = [];
for jt = 1:length(chck)
    selIdx(jt) = ~isempty(chck{jt});
end
selIdx = find(selIdx~=0);
scalef = str2double(hdrCSC{selIdx}(min(regexp(hdrCSC{selIdx},'\d')):end));


dataSamplesCSC = double(dataSamplesCSC(:))'; % raw data in samples
dataSamplesCSC = dataSamplesCSC.*scalef.*1e6;

% setting up the filter (300 to 3000Hz)
[a, b] = butter(5, [0.01875 0.1875], 'bandpass');
filtSig = filtfilt(a,b,dataSamplesCSC);

% start and end times
startSig = 32000*1000;
endSig = size(dataSamplesCSC,2);

% visualize raw data
figure(1)
handle1 = plot(dataSamplesCSC(startSig:endSig),'color',[0 0 1]);
% xticks([0:1:120])
% xlabel('Time in Seconds');
ylabel('\muV', 'FontSize', 16);
% set(gca,'YTickLabel',[]);
set(gca, 'YLim', [-750, 750]);
set(gca, 'XLim', [0, endSig-startSig], 'XTick', 0:3200:endSig-startSig,...
    'XTickLabel', 0:0.1:1);
mhandle = gca;
mhandle.FontSize = 20;
% set(gca,'xtick',[])
set(gca,'visible','off')


% visualize filtered data
figure(2)
handle2 = plot(filtSig(startSig:endSig),'color',[0.75 0.75 0.75]);
xlabel('Time in Seconds');
ylabel('\muV', 'FontSize', 16);
% set(gca,'YTickLabel',[]);
set(gca, 'YLim', [-80, 80]);
set(gca, 'XLim', [0, endSig-startSig], 'XTick', 0:3200:endSig-startSig,...
    'XTickLabel', 0:0.1:1);
mhandle = gca;
mhandle.FontSize = 20;
set(gca,'visible','off')

% LFP
[a2, b2] = butter(1, [0.00003125 0.0187], 'bandpass');
filtSig_lfp = filtfilt(a2,b2,dataSamplesCSC);
figure(3)
plot(filtSig_lfp(startSig:endSig),'color',[0 0 0]);
ylabel('\muV', 'FontSize', 16);
set(gca, 'YLim', [-750, 750]);
set(gca, 'XLim', [0, endSig-startSig], 'XTick', 0:3200:endSig-startSig,...
    'XTickLabel', 0:0.1:1);
mhandle = gca;
mhandle.FontSize = 20;
set(gca,'visible','off')

% [0.105882353 0.678431373 0.811764706];
%% RAW WAVESHAPE
startP = 1009*32000;
endP = 1010*32000;
cl = 1;

% specific cluster (cl)
idx = cluster_class(:,1)==cl & cluster_class(:,2)>=startP & cluster_class(:,2)<=endP; % 1588

% all cluster
idx = cluster_class(:,2)>=startP & cluster_class(:,2)<=endP;

relClust = spikes(idx,:);

% visualize cluster "cl" within time window
figure(1)
hold on
for ix = 1:size(relClust,1)
    plot(relClust(ix,:),'r');
end
% hold off

idx = cluster_class(:,1)==cl & cluster_class(:,2)<=startP | cluster_class(:,1)==cl & cluster_class(:,2)>=endP; % 276
othClust = spikes(idx,:);

% visualize cluster "cl" outside time window
% figure(2)
% hold on
for ix = 1:size(othClust,1)
    plot(othClust(ix,:),'r');
end
hold off
