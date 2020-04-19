%% tEMt bham pipeline
% load dependencies
addpath(genpath('X:\hanslmas-ieeg-compute\Luca\functions\wave_clus-master')); % waveclus 3.0
addpath('X:\hanslmas-ieeg-compute\Luca\functions'); % my functions
addpath('X:\hanslmas-ieeg-compute\Luca\functions\Neuralynx_19012019'); % Neuralynx (the commons folder function doesnt work on my PC)
% addpath(genpath('X:\hanslmas-ieeg-compute\Luca\TREBER\Scripts'))

%% Start tuningA - PRE
cd('C:\Experiments\tuning_tEMt_041219\')
edit startConceptTuning.m
% change/check the required parameter (prePost = 'pre' ;session; subject; trg.debug for a practice run, else:
% trg.ttl)

%% Spikedetection and Automatic spikesorting
% copy the macro .ncs and logfiles here
calcFolder = 'C:\Experiments\tuning_tEMt_041219\sub-0014\dat_sesh1\'; % dont forget the backslash here!!

mkdir(calcFolder);
cd(calcFolder);

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
Do_clustering('all', 'parallel', true, 'make_plots', false);
% Do_clustering('CSC_antHippR1_spikes.mat', 'parallel', false,'make_plots', true);

% % negative clustering
cd ..\negDetect
Do_clustering('all', 'parallel', true, 'make_plots', false);

%% load all sorted cluster into variable "allCl"
% load negative cluster
allCl = {};
negCl = dir('times_*.mat');
for cl = 1:size(negCl,1)
    load(negCl(cl).name, 'cluster_class')
    test = negCl(cl).name
    test = regexprep(test,'times_CSC_', '');
    test(end-3:end ) = [];
    test = strjoin(test,'neg', num2str(cl));
    for wc = 1:max(cluster_class(:,1))
        clustTimes = cluster_class(cluster_class(:,1)==wc ,2)./32000;
        allCl = [allCl {clustTimes}];
    end
end

% load positive cluster
cd ..\posDetect\
posCl = dir('times_*.mat');
for cl = 1:size(posCl,1)
    load(posCl(cl).name, 'cluster_class')
    for wc = 1:max(cluster_class(:,1))
        clustTimes = cluster_class(cluster_class(:,1)==wc ,2)./32000;
        allCl = [allCl {clustTimes}];
    end
end
  
%% load logfile
cd(calcFolder)
[mCell, allTrials, trlTrigger] = loadLogs_tEMt(calcFolder);
locking = 1;

% creates a SUxTRLs cell with the spikestamps of each SU in each trl
timeWindow = [0.3 1];
allSU = {};
for cl = 1:size(allCl,2)
    spikeTimesCluster = allCl{cl};
    spiketimes_clusterTrial = insertSpiketimes(trigger, spikeTimesCluster, locking, timeWindow);
    allSU = [allSU; spiketimes_clusterTrial];
end

% same as above, but for the baseline period
allBL = {};
timeWindow = [-1 -0.3];
for cl = 1:size(allCl,2)
    spikeTimesCluster = allCl{cl};
    spiketimes_clusterTrial = insertSpiketimes(trigger, spikeTimesCluster, locking, timeWindow);
    allBL = [allBL; spiketimes_clusterTrial];
end

% transform spiketimes to spikenumber
allBLumber = cellfun(@length,allBL);
meanBL = mean(allBLumber,2);
stdBL = std(allBLumber,0,2);
allSUmber = cellfun(@length,allSU);

% threshold for concept activity
th = meanBL+5*stdBL;

%% collapse over repetitions of the same image and compare against th
concImg = {};
for ix = 1:size(mCell,2)
    trls = mCell{2,ix}; % those are the trials in which the same image is shown
    empDat = median(allSUmber(:,trls),2); % calculate the median firing rate towards that image over all repetition for each SU. Results in a SUx1 vector
    
    if any(empDat>th & empDat>=2) % if any of the single unit is a concept neuron to the image, kill the image
        concImg = [concImg, mCell(1,ix)];
    end
end

%% copy the images that were not identified in the tuning session to tEMt to be run for the EM part
switch sesh
    case '01'
        sourceFolder = 'C:\Experiments\tuning_tEMt_041219\EMpairs_v5_2017-05-01\image_data\Tune\sesh1';
    case '02'
        sourceFolder = 'C:\Experiments\tuning_tEMt_041219\EMpairs_v5_2017-05-01\image_data\Tune\sesh2';
end
seshFolder = 'C:\Experiments\tEMt_041219\EMpairs_v5_2017-09-11\image_data\EMtune\formatted_180-180';
cd(sourceFolder)
ImgFold = dir('*.jpg');

killpic = ismember({ImgFold.name},concImg);

for iy = 1:size(ImgFold,1)
    if killpic(iy) == 1
        continue
    end
    copyfile( [sourceFolder,filesep,ImgFold(iy).name], seshFolder);
end

%% start EM task
restoredefaultpath
cd('C:\Experiments\tEMt_041219');
edit startEMexp.m
% run the experiment from the startEMexp script after youve checked the parameters

%% start TuningB - POST
restoredefaultpath
cd('C:\Experiments\tuning_tEMt_041219\')
edit startConceptTuning.m
% change/check the required parameter (prePost = 'pre' ;session; subject; trg.debug for a practice run, else:
% trg.ttl)

