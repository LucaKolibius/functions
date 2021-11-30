% housekeeping: George mentioned that microrecordings are clipped. sub-1003
% is sampled at 1024. Maybe redo MW (including interpolation and AR) myself

% spike-ripple coincidences in IU during indexing vs. during normal trials vs. coincidences in SU/GU inside indexing bundle vs coincidences in SU/GU (outside indexing bundle)
% hits vs misses (currently only hits) in encoding and retrieval
% don't average over SU!
% differentiate between indexed trials in gray units
% ripple power / frequency / duration / timing?
% only ripples from MW instead of bundle

closeVars; clear; clc;
load('X:\Luca\data\allSbj\rpplRec.mat', 'allSU');
load('X:\Luca\data\allSbj\allTrig.mat', 'trigALL', 'trigIU');

% 113 bundles

%% compare the number of spike - ripple conincidence between index units and normal single units (excluding gray untis)
coincSU = [];
coincGU = [];
coincIU = [];

eep = [];
for su = 1 : size(allSU,2)
    
    % IF WE HAVE NO LFP / RIPPLE INFO, SKIP
    if isempty(allSU(su).su)
        continue
    end
    
    rpplPerc                = mean(allSU(su).rpplRec);                                            % the percantage of the whole recording that is covered in ripples
    spksDurRppl             = mean(allSU(su).rpplRec(allSU(su).spks(:)) ); % the percentage of all spikes that occur during ripples
   
    
    % 2999561
    % 2934110
    % su = 446
    
    switch trigALL(su).iu
        case 0
            coincSU             = [coincSU, spksDurRppl / rpplPerc];  % a value >1 indicates that this IU spikes more often during ripples
        case 1
            coincGU             = [coincGU, spksDurRppl / rpplPerc];
        case 2
            coincIU             = [coincIU, spksDurRppl / rpplPerc];
    end
end

dt = 0:0.5:22;
close all
figure('units','normalized','outerposition',[0 0 1 1])

subplot(311); hold on
histogram(coincIU, dt, 'Normalization','probability');
plot([nanmedian(coincIU) nanmedian(coincIU)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');

subplot(312); hold on
histogram(coincGU, dt, 'Normalization','probability');
plot([nanmedian(coincGU) nanmedian(coincGU)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');

subplot(313); hold on
histogram(coincSU, dt, 'Normalization','probability');
plot([nanmedian(coincSU) nanmedian(coincSU)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');

%% RIPPLE QUANTITY AND COINCIDENCES
%  split up into IU (idxd trl or not) / GU (idxd trl or not) / SU (on IU bundle or not)
closeVars; clearvars -except allSU trigALL trigIU; clc;

if ~exist('allSU','var'); load('X:\Luca\ripple and IU\ripplRec.mat', 'allSU'); end
load('X:\Luca\data\allSbj\allTrig.mat', 'trigALL', 'trigIU');
tw = 0; % 0 = whole trial | 1 = preCue | 2 = periCue | 3 = periResp

[resQuant, resCoinc, resTime] = spkRppl_anal_sub1 (allSU, trigALL, trigIU, tw);

% VISUALIZE
sets.tw      = tw;
sets.coinc   = 0; % 0 = QUANTITY | 1 = COINCIDENCE
sets.avgTrl  = 0; % 0 for each trial independently | 1 for average over one SU
sets.del0    = 0; % 0 = keep 0s  | 1 = delete 0s (this isn't implemented for avgSU - either kick out the trials with 0 or the SU with 0...)
sets.medmean = 1; % 0 for mean | 1 for median

spkRppl_visu(resQuant, resCoinc, resTime, sets)
saveas(gca, 'visu.png')

% spike-ripple coincidences in IU during indexing vs. during normal trials vs. coincidences in SU/GU inside indexing bundle vs coincidences in SU/GU (outside indexing bundle)  
% hits vs misses (currently only hits) in encoding and retrieval 
% don't average over SU!
% differentiate between indexed trials in gray units
% ripple power / frequency / duration / timing?
% only ripples from MW instead of bundle
