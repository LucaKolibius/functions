for i = 1 : size(allSubj,1)
    subjID = allSubj{i};
    SME(subjID)
end

[gausConvENC_hitsALL, gausConvENC_missALL, ~, ~] = extractSME(allSubj); % SME for each session

% mean SME over all sessions
gausConvENC_hitsALL_m = nanmean(gausConvENC_hitsALL,1);
gausConvENC_missALL_m = nanmean(gausConvENC_missALL,1);

% standard error = std/sqrt(n)
gausConvENC_hitsALL_se = std(gausConvENC_hitsALL,1) / sqrt(size(allSubj,1));
gausConvENC_missALL_se = nanstd(gausConvENC_missALL,1) / sqrt(size(allSubj,1));

% visualize
visuSME

%% inference statistics
% the kernel is already taken away / steps of 0.0002
% 1:4625; % fix to cue; 1s - kernel
% 4626:14625; % cue to stim; 2s
% 14626:29250; % stim to resp - kernel; 3s

% hits
encHits_fixCue   = gausConvENC_hitsALL(: , 1:4625); % -1 : 0
encHits_cueStim  = gausConvENC_hitsALL(: , 4626:14625); % 0 : 2
encHits_stim1    = gausConvENC_hitsALL(: , 14626:19625); % 2 : 3
encHits_stimResp = gausConvENC_hitsALL(: , 19626:29250); % 3 : 5

% misses
encMiss_fixCue   = gausConvENC_missALL(:, 1:4625);
encMiss_cueStim  = gausConvENC_missALL(:, 4626:14625);
encMiss_stim1    = gausConvENC_missALL(:, 14626:19625);
encMiss_stimResp = gausConvENC_missALL(:, 19626:29250);

% means: hits
encHits_fixCue_m   = nanmean(encHits_fixCue, 2);
encHits_cueStim_m  = nanmean(encHits_cueStim, 2);
encHits_stim1_m    = nanmean(encHits_stim1, 2);
encHits_stimResp_m = nanmean(encHits_stimResp, 2);

% means: miss
encMiss_fixCue_m   = nanmean(encMiss_fixCue, 2);
encMiss_cueStim_m  = nanmean(encMiss_cueStim, 2);
encMiss_stim1_m    = nanmean(encMiss_stim1, 2);
encMiss_stimResp_m = nanmean(encMiss_stimResp, 2);

% hits minus misses (average)
encDiff_fixCue_m   = encHits_fixCue_m   - encMiss_fixCue_m;
encDiff_cueStim_m  = encHits_cueStim_m  - encMiss_cueStim_m;
encDiff_stim1_m    = encHits_stim1_m    - encMiss_stim1_m;
encDiff_stimResp_m = encHits_stimResp_m - encMiss_stimResp_m;

forSPSShimi = [];
forSPSSdiff = [];
forSPSShimi = [encHits_fixCue_m, encHits_cueStim_m, encHits_stim1_m, encHits_stimResp_m, ones(size(encHits_fixCue_m,1),1); encMiss_fixCue_m, encMiss_cueStim_m, encMiss_stim1_m, encMiss_stimResp_m, zeros(size(encMiss_fixCue_m,1),1)];
forSPSSdiff = [encDiff_fixCue_m, encDiff_cueStim_m, encDiff_stim1_m, encDiff_stimResp_m];

clearvars -except forSPSShimi forSPSSdiff
open forSPSShimi
open forSPSSdiff