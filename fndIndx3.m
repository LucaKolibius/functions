% ideas:
% segmentierung von spike detection sollte auf varianz basieren!
addpath(genpath('X:\hanslmas-ieeg-compute\Luca_old\functions\wave_clus-master')); % waveclus 3.0
addpath('X:\hanslmas-ieeg-compute\Luca_old\functions'); % my functions
addpath('X:\hanslmas-ieeg-compute\Luca_old\functions\Neuralynx_19012019'); % Neuralynx (the commons folder function doesnt work on my PC)
addpath(genpath('X:\Luca_old\TREBER\Scripts'))


try
p2d = cd;
p2d(end+1)='\';
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, hitsIdx, missIdx, allTrials, sessionNum, retTrigger, encTrigger]=loadLogs(p2d);
dts.start = encTrigger(1);
dts.end = max(max(retTrigger))+1; % remember that retTrigger is not linearly increasing
dts.dt = 0.001:0.001:dts.end; % create a linerally spaced vector from the beginning to the end of the experiment in steps of 1ms

catch
    disp('Using default dt values');
    dts.start = 1;
    dts.end = 4000;
    dts.dt = 3;
end

% positive detection
par = set_parameters();
par.detection = 'pos';
Get_spikes('all', 'parallel', true, 'par', par);
movefile('*.mat', 'posDetect'); % moves all .mat files that include the just detected spikes into the folder "posDetect"

% % negative detection
par.detection = 'neg';
Get_spikes('all', 'parallel', true, 'par', par);
movefile('*.mat', 'negDetect'); % moves all .mat files that include the just detected spikes into the folder "negDetect"

% positive clustering
cd posDetect
Do_clustering('all', 'parallel', true, 'make_plots', true);
% Do_clustering('CSC_antHippR1_spikes.mat', 'parallel', false,'make_plots', true);

% % negative clustering
cd ..\negDetect
Do_clustering('all', 'parallel', true, 'make_plots', true);

clear it
% positive manual inspection
cd ..\posDetect
if ~exist('it'); it=1; end
manuSorting(dts, it);
clear it

%     checkchannelLDK(it)
%     mkSpiketimes(dts ,it);

%     pnClus = dir('times_*.mat');
%     wave_clus(pnClus(it).name);
%     load('it.mat'); disp(it);

% negative manual inspection
clear it; cd ..\negDetect;
if ~exist('it'); it=1; end
manuSorting(dts, it);
clear it


%%
clear all
cd('X:\Luca\sub-1013');
allSesh = dir('2020*');

for ia = 1:length(allSesh)
    cd([allSesh(ia).folder, filesep allSesh(ia).name])
    
%     % positive detection
%     par = set_parameters();
%     par.detection = 'pos';
%     Get_spikes('all', 'parallel', true, 'par', par);
% %     Get_spikes('Micro_MRTC8_.ncs', 'parallel', false, 'par', par);
%     movefile('*.mat', 'posDetect'); % moves all .mat files that include the just detected spikes into the folder "posDetect"
%     
%     % % negative detection
%     par.detection = 'neg';
%     Get_spikes('all', 'parallel', true, 'par', par);
%     movefile('*.mat', 'negDetect'); % moves all .mat files that include the just detected spikes into the folder "negDetect"
    
    % positive clustering
    cd posDetect
    Do_clustering('all', 'parallel', false, 'make_plots', true);
    % Do_clustering('CSC_antHippR1_spikes.mat', 'parallel', false,'make_plots', true);
    
    % % negative clustering
    cd ..\negDetect
    Do_clustering('all', 'parallel', false, 'make_plots', true);
end


%% load lfp
%% define the FieldSelection variable and set the ExtracMode.
FieldSelection(1) = 1;%timestamps
FieldSelection(2) = 0;
FieldSelection(3) = 0;%sample freq
FieldSelection(4) = 0;
FieldSelection(5) = 1;%samples
ExtractHeader = 1;
ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.

[p2d] = 'X:/';
cd('X:\Luca\toolboxes\Neuralynx_19012019\') % after restarting the neuralynx function does not find the path to a mex file, unless you cd into that folder

[CSCfiles] = dir([p2d,'*.ncs']);

dum = cell(1,length(CSCfiles));
it = 1;

[timestampsCSC, dataSamplesCSC,hdrCSC] = Nlx2MatCSC([p2d,CSCfiles(it).name], FieldSelection, ExtractHeader, ExtractMode, []);% import the raw data

%extract the scale factor
chck   = regexp(hdrCSC,'ADBitVolts');
selIdx = [];
for jt = 1:length(chck); selIdx(jt) = ~isempty(chck{jt});end;
selIdx = find(selIdx~=0);
scalef = str2double(hdrCSC{selIdx}(min(regexp(hdrCSC{selIdx},'\d')):end));

% extract the sampling-frequency
chck = regexp(hdrCSC,'SamplingFrequency');
selIdx = [];
for jt = 1:length(chck); selIdx(jt) = ~isempty(chck{jt});end;
selIdx = find(selIdx~=0);
Fs = str2double(hdrCSC{selIdx}(min(regexp(hdrCSC{selIdx},'\d')):end));

%flatten
dataSamplesCSC=double(dataSamplesCSC(:))';
dataSamplesCSC = dataSamplesCSC.*scalef.*1e6;
% create a time vector of theCSCtime
[CSCTime] = 0:1/Fs:(length( dataSamplesCSC )-1)/Fs; % in units of sec
dum{it} = [];
dum{it}.trial{1} = dataSamplesCSC;
dum{it}.time{1} = CSCTime;
dum{it}.label = {chanLab};
dum{it}.cfg.hdr = hdrCSC;
dum{it}.fsample = Fs;

mTime = size(CSCTime,2); 
mTime = mTime/32