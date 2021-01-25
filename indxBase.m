% % This should be my new base script
% fetch_allSubj;

whereAmI(1) % own laptop

% set paths
addpath(genpath('X:\Luca\functions\wave_clus-master')); % waveclus 3.0
% addpath('X:\Luca\functions\Neuralynx_19012019'); % Neuralynx (the commons folder function doesnt work on my PC)
addpath(genpath('X:\Luca\TREBER\Scripts'))
addpath('X:\Common\fieldtrip-20200310');
addpath('X:\Common\toolboxes\fieldtrip-20200310');
ft_defaults;
dbstop if error
dbclear if error

addpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\toolboxes\fieldtrip-20200603')
addpath(genpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\functions')); % my functions
fixDoubleLFP % fix: on the same bundle IU1 indexes trl1 and IU2 does not index trl1. don't count lfp double.
findIndexSU_hz(99, 1.645, 1.645, 0, 0) % find index units
lfp2microLFP % save LFP from georges structure into mine 
lfp2rpplBund
rppl_statsVisu
preCuePowDiff % precue Power || try out uniqFreq per freqBand & try out time-frequency
preCuePowDiff_perm % statistics on it
spkRppl_ppc
spkPhase_encRet % favChan

TW = [1000 5000]; % 1s before and 5s after
trigCode = [1]; % cue-locked
TFtoi3(trigCode, TW); %% TIME FREQUENCY ANALYSIS

TW = [5000 1000]; % 5s before and 1s after
trigCode = [3]; % resp-locked
TFtoi3(trigCode, TW);

TW = [1000 1000]; % 1s before and 1s after
trigCode = [1 3]; % cue- & resp-locked
TFtoi3(trigCode, TW);

%%
% compare micro LFP of precue period (-1:0) for indexed and non-indexed
% trials
LFPanal

%% detect and cluster spikes (manually and automatically) using fndIndx3.m

%  Here I use cross-correlation to get rid of some artefacts
subjID = 'P11ERL_S2';
subjID = 'P0#_S#';
disp(subjID);

xc_postClust(subjID); % does the crosscorrelation
quickfix = xc_postClust2(subjID); % rejects cluster based on crosscorrelation
if ~isempty(quickfix); open quickfix; end
rename_clNames(subjID); % fixes possible label errors

[encSpiketimes_cueLocked, encSpiketimes_respLocked, encSpiketimes_preTrial, encSpiketimes_stimLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); % put spiketimes into the table template, output to visually check
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, encSpiketimes_preTrial, encSpiketimes_stimLocked, clNames, mCD, subjID);

% extract results for each participant in a way that works for SPSS
extractMeanRSA

% generate a mean RSA for all comparisons for hits as well as for misses that can be visualized
Average4RSA(allSubj)

% visualize all single units
visuSU_old % visualizes hits and misses for encoding and retrieval + frequency plot + waveshape
visuSU('X:\Luca\indexSUvisu\varWindow_ENC_0-213_RET_-1113_DP_p99_th1')

% compares within trial (main diagnonal) and within category (cross diagonal) and outputs values as SPSS input
RSA_crossidag_WtWc % currently for hippocampus cuecue hit+miss

% SME
% script that runs function "SME" for all subjects (allSubj)
% then it runs "extractSME". This function summarizes the output from "SME" for all subjects (allSubj)
runSME

% visualizes WT and WC (based on crossdiagnonal) as a bar plot with within
% subject lines
% resides in luca\saved functions
% addpath('X:\Luca\saved functions')
WIP_visu_250919

% finds significantly high within-trial values of the RSA
fndSignWT

% based on permutation testing determines significant correlations between
% the TS of a single unit in encoding and retrieval. Then based on
% permutation testing identifies the number of expected significant
% correlations based on the H0
encRetTSall

% same thing but without the correlation (only DP)
findIndexSU(99, 1, 0.524)
findIndexSU_varWind(99,1)

% looks for a significant autocorrelation for a single unit for the TS of
% encoding or retrieval. Optional to use retrieval at presentation order
tempGrad_AuC_encRet % temporal gradient autocorrelation encoding retrieval

% desktop
% looking at XC of spiketimes between encoding and retrieval for indexed
% trials
xcSpikes_encret

% analysis of tEMt data
anal_sub_1013

% find and extract LFP for index neurons (done)
% for ic = [4 6:13] % takes overlap into account
fndLFP(ic)
% end

% converts the timestamps from the fVSp structure into the allSU strucutre
% considering electrode location (elecLoc) & if its a good unit (clNames)
conv2allSU

%% script I use to find out whether index neurons fire first during encoding or during retrieval
indexNeuron_timeLag_retEnc

%% to do
slyWindow % not finished

% permutation test on number of significant windows
WIP_mar2020

%%
clear; clc
fetch_allSubj;
for sbj = 1 : size(allSubj, 1)
subjID = allSubj{sbj};
disp(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, encSpiketimes_preTrial, encSpiketimes_stimLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); % put spiketimes into the table template, output to visually check
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, encSpiketimes_preTrial, encSpiketimes_stimLocked, clNames, mCD, subjID);
end


%%
WIP_visu_250919 % visualisation of WT and BT for hippocampus CueCue hits and misses
fndSignWT % finds significant correlations in the hit within trial (WT) RSA for CueCue in the hippocampus
fndSignBT % same for between trial (BT)
encRetTSall % all SU time series between encoding and retrieval


