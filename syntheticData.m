clear

% possible stimulus combinations
stimu{1} = 'pp';
stimu{2} = 'ff';
stimu{3} = 'fp';

% is it a hit or a miss?
hitMiss{1} = 'hit';
hitMiss{2} = 'miss';

% number of trials
allTrials = 200;

% for each session
for sesh = 1:40
    clearvars -except stimu hitMiss allTrials sesh
    subjID = strcat('P99_S', num2str(sesh));

    % randomly select stimulus combination and hit/miss for each trial
for trl = 1:allTrials
    trial(1,trl)     = stimu(randi(3,1));
    trial(2,trl)     = hitMiss(randi(2,1));  
end


% create hit and miss index for later analysis
hitsIdx = zeros(allTrials,1);
missIdx = zeros(allTrials,1);
for trl = 1:allTrials
    if strcmp(trial(2,trl),'hit')
        hitsIdx(trl) = 1;
        missIdx(trl) = 0;
    else
        hitsIdx(trl) = 0;
        missIdx(trl) = 1;
    end
end

temp = contains(trial(2,:),'miss');
missIdx  = find(temp==1)';

temp = contains(trial(2,:),'hit');
hitsIdx = find(temp==1)';

% spike number per SU
trialEnc = trial;
trialRet = trial;
ff = [5 10 20 30]; % possible average spike numbers
sigm = ff./5; % variance

% we have 30 SU
for su = 1:30
    rndidx = randi(4); % random index
    x = -1;
    
    while any(x < 0) % this is the spike number vector for one SU, make sure we do not have negative numbers, by resampling
        x = round( normrnd( ff(rndidx), sigm(rndidx), 1, allTrials*2));
    end
    
    trialEnc(su+2, 1:allTrials) = num2cell(x(1:allTrials)); % the first 200 columns belong to encoding
    trialRet(su+2, 1:allTrials) = num2cell(x(allTrials+1:allTrials*2)); % the last 200 coluns belong to retrieval
end

% make into a table
encSpikeNumber_syn = cell2table(trialEnc);
retSpikeNumber_syn = cell2table(trialRet);

%%
try
    cd Z:/Luca/data
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
mkdir(mSession)
cd(mSession)
mkdir('whatever');
cd whatever
mkdir('advancedAnalysis');
cd advancedAnalysis\


%% number of FF/PP/FP for hits & misses
% hit
numFF_hit=0;
for i=1:size(encSpikeNumber_syn,2)
    if strcmp(encSpikeNumber_syn{1,i}, 'ff') && strcmp(encSpikeNumber_syn{2,i}, 'hit')
        numFF_hit=numFF_hit+1;
    end
end

numPP_hit=0;
for i=1:size(encSpikeNumber_syn,2)
    if strcmp(encSpikeNumber_syn{1,i}, 'pp') && strcmp(encSpikeNumber_syn{2,i}, 'hit')
        numPP_hit=numPP_hit+1;
    end
end
numPP_hit=numPP_hit+numFF_hit;

numFP_hit=0;
for i=1:size(encSpikeNumber_syn,2)
    if strcmp(encSpikeNumber_syn{1,i}, 'fp') && strcmp(encSpikeNumber_syn{2,i}, 'hit')
        numFP_hit=numFP_hit+1;
    end
end
numFP_hit=numFP_hit+numPP_hit;
numHit=[numFF_hit numPP_hit numFP_hit];

% miss
numFF_miss=0;
for i=1:size(encSpikeNumber_syn,2)
    if strcmp(encSpikeNumber_syn{1,i}, 'ff') && strcmp(encSpikeNumber_syn{2,i}, 'miss')
        numFF_miss=numFF_miss+1;
    end
end

numPP_miss=0;
for i=1:size(encSpikeNumber_syn,2)
    if strcmp(encSpikeNumber_syn{1,i}, 'pp') && strcmp(encSpikeNumber_syn{2,i}, 'miss')
        numPP_miss=numPP_miss+1;
    end
end
numPP_miss = numPP_miss+numFF_miss;

numFP_miss = 0;
for i=1:size(encSpikeNumber_syn,2)
    if strcmp(encSpikeNumber_syn{1,i}, 'fp') && strcmp(encSpikeNumber_syn{2,i}, 'miss')
        numFP_miss=numFP_miss+1;
    end
end
numFP_miss = numFP_miss+numPP_miss;
numMiss = [numFF_miss numPP_miss numFP_miss];

% quickly save the number of trials, hits and misses
mHits = numHit(3);
mMiss = numMiss(3);
hitsMiss = [mHits+mMiss mHits mMiss];
save(['hitsMiss', subjID, '.mat'], 'hitsMiss');

%
enoughC    = 1;
enoughO    = 1;
enoughH    = 1;
enoughMiss = 1;

%% extract spiketimes
% cortex
encSpikeNumber_cueLocked_c   = encSpikeNumber_syn;
encSpikeNumber_respLocked_c  = encSpikeNumber_syn;
encSpikeNumber_preTrial_c    = encSpikeNumber_syn;
encSpikeNumber_stimLocked_c  = encSpikeNumber_syn;
retSpikeNumber_cueLocked_c   = retSpikeNumber_syn;
retSpikeNumber_respLocked_c  = retSpikeNumber_syn;

% hipp
encSpikeNumber_cueLocked_h   = encSpikeNumber_syn;
encSpikeNumber_respLocked_h  = encSpikeNumber_syn;
encSpikeNumber_preTrial_h    = encSpikeNumber_syn;
encSpikeNumber_stimLocked_h  = encSpikeNumber_syn;
retSpikeNumber_cueLocked_h   = retSpikeNumber_syn;
retSpikeNumber_respLocked_h  = retSpikeNumber_syn;

% other
encSpikeNumber_cueLocked_o   = encSpikeNumber_syn;
encSpikeNumber_respLocked_o  = encSpikeNumber_syn;
encSpikeNumber_preTrial_o    = encSpikeNumber_syn;
encSpikeNumber_stimLocked_o  = encSpikeNumber_syn;
retSpikeNumber_cueLocked_o   = retSpikeNumber_syn;
retSpikeNumber_respLocked_o  = retSpikeNumber_syn;

mkdir('elecLoc'); cd elecLoc;
mkdir('spikenumbers'); cd spikenumbers;
save(['spikeNumber_hipp_', subjID, '.mat'],  'encSpikeNumber_cueLocked_h', 'encSpikeNumber_respLocked_h', 'encSpikeNumber_preTrial_h', 'encSpikeNumber_stimLocked_h', 'retSpikeNumber_cueLocked_h', 'retSpikeNumber_respLocked_h');
save(['spikeNumber_cort_', subjID, '.mat'],  'encSpikeNumber_cueLocked_c', 'encSpikeNumber_respLocked_c', 'encSpikeNumber_preTrial_c', 'encSpikeNumber_stimLocked_c', 'retSpikeNumber_cueLocked_c', 'retSpikeNumber_respLocked_c');
save(['spikeNumber_other_', subjID, '.mat'], 'encSpikeNumber_cueLocked_o', 'encSpikeNumber_respLocked_o', 'encSpikeNumber_preTrial_o', 'encSpikeNumber_stimLocked_o', 'retSpikeNumber_cueLocked_o', 'retSpikeNumber_respLocked_o');

%% hipp
cd ..\..
mkdir('RSA');
cd RSA
mkdir('oldSize');
cd oldSize
mkdir('hippocampus');
mkdir('cortex')
mkdir('other');
cd hippocampus

if enoughH == 1
    %% Hits (Cue/Cue) - 1
    [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_h, retSpikeNumber_cueLocked_h, 'hit', 0);
    RSAhits_CueCue_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h); % normalizes spikenumber and correlates enc+ret for RSA
    
    % Misses (Cue / Cue)
    if enoughMiss == 1 % only calculate this if there are enough misses
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_h, retSpikeNumber_cueLocked_h, 'miss', 0);
        RSAmiss_CueCue_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    else
        RSAmiss_CueCue_h = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Cue / Cue)';
    plotRSA(RSAhits_CueCue_h, RSAmiss_CueCue_h, comp, hitsIdx, missIdx);
    savefig(['RSA_hipp_CueCue', subjID]);
    saveas(gcf,['RSA_hipp_CueCue' subjID,'.png']);
    
    %% Hits (Cue/Resp) - 2
    [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_h, retSpikeNumber_respLocked_h, 'hit', 0);
    RSAhits_CueResp_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    
    % CueResp Misses
    if enoughMiss == 1
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_h, retSpikeNumber_respLocked_h, 'miss', 0);
        RSAmiss_CueResp_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    else
        RSAmiss_CueResp_h = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Cue / Resp)';
    plotRSA(RSAhits_CueResp_h, RSAmiss_CueResp_h, comp, hitsIdx, missIdx);
    savefig(['RSA_hipp_CueResp', subjID]);
    saveas(gcf,['RSA_hipp_CueResp' subjID,'.png']);
    
    %% Hits (Resp/Cue) - 3
    [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_respLocked_h, retSpikeNumber_cueLocked_h, 'hit', 0);
    RSAhits_RespCue_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    
    % RespCue Misses
    if enoughMiss == 1
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_respLocked_h, retSpikeNumber_cueLocked_h, 'miss', 0);
        RSAmiss_RespCue_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    else
        RSAmiss_RespCue_h = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Resp / Cue)';
    plotRSA(RSAhits_RespCue_h, RSAmiss_RespCue_h, comp, hitsIdx, missIdx);
    savefig(['RSA_hipp_RespCue', subjID]);
    saveas(gcf,['RSA_hipp_RespCue' subjID,'.png']);
    
    %% Hits (Resp / Resp) - 4
    [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_respLocked_h, retSpikeNumber_respLocked_h, 'hit', 0);
    RSAhits_RespResp_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    
    % RespResp Misses
    if enoughMiss == 1
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_respLocked_h, retSpikeNumber_respLocked_h, 'miss', 0);
        RSAmiss_RespResp_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    else
        RSAmiss_RespResp_h = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Resp / Resp)';
    plotRSA(RSAhits_RespResp_h, RSAmiss_RespResp_h, comp, hitsIdx, missIdx);
    savefig(['RSA_hipp_RespResp', subjID]);
    saveas(gcf,['RSA_hipp_RespResp' subjID,'.png']);
    
    %% Hits (Pre / Cue) - 5
    [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_preTrial_h, retSpikeNumber_cueLocked_h, 'hit', 0);
    RSAhits_PreCue_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    
    % PreCue Misses
    if enoughMiss == 1
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_preTrial_h, retSpikeNumber_cueLocked_h, 'miss', 0);
        RSAmiss_PreCue_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    else
        RSAmiss_PreCue_h = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Pre / Cue)';
    plotRSA(RSAhits_PreCue_h, RSAmiss_PreCue_h, comp, hitsIdx, missIdx);
    savefig(['RSA_hipp_PreCue', subjID]);
    saveas(gcf,['RSA_hipp_PreCue' subjID,'.png']);
    
    %% Hits (Pre / Resp) - 6
    [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_preTrial_h, retSpikeNumber_respLocked_h, 'hit', 0);
    RSAhits_PreResp_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    
    % PreResp Misses
    if enoughMiss == 1
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_preTrial_h, retSpikeNumber_respLocked_h, 'miss', 0);
        RSAmiss_PreResp_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    else
        RSAmiss_PreResp_h = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Pre / Resp)';
    plotRSA(RSAhits_PreResp_h, RSAmiss_PreResp_h, comp, hitsIdx, missIdx);
    savefig(['RSA_hipp_PreResp', subjID]);
    saveas(gcf,['RSA_hipp_PreResp' subjID,'.png']);
    
    %% Hits (Stim / Cue) - 7
    [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_h, retSpikeNumber_cueLocked_h, 'hit', 0);
    RSAhits_StimCue_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    
    % StimCue Misses
    if enoughMiss == 1
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_h, retSpikeNumber_cueLocked_h, 'miss', 0);
        RSAmiss_StimCue_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    else
        RSAmiss_StimCue_h = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Stim / Cue)';
    plotRSA(RSAhits_StimCue_h, RSAmiss_StimCue_h, comp, hitsIdx, missIdx);
    savefig(['RSA_hipp_StimCue', subjID]);
    saveas(gcf,['RSA_hipp_StimCue' subjID,'.png']);
    
    %% Hits (Stim / Resp) - 8
    [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_h, retSpikeNumber_respLocked_h, 'hit', 0);
    RSAhits_StimResp_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    
    % StimResp Misses
    if enoughMiss == 1
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_h, retSpikeNumber_respLocked_h, 'miss', 0);
        RSAmiss_StimResp_h = normSpikeNumber(tempCell_enc_h, tempCell_ret_h);
    else
        RSAmiss_StimResp_h = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Stim / Resp)';
    plotRSA(RSAhits_StimResp_h, RSAmiss_StimResp_h, comp, hitsIdx, missIdx);
    savefig(['RSA_hipp_StimResp', subjID]);
    saveas(gcf,['RSA_hipp_StimResp' subjID,'.png']);
    
    %% save all RSAs
    save(['RSAmat_hipp', subjID, '.mat'], 'RSAhits_CueCue_h', 'RSAmiss_CueCue_h', 'RSAhits_CueResp_h', 'RSAmiss_CueResp_h', 'RSAhits_RespCue_h', 'RSAmiss_RespCue_h', 'RSAhits_RespResp_h', 'RSAmiss_RespResp_h', 'RSAhits_PreCue_h', 'RSAmiss_PreCue_h', 'RSAhits_PreResp_h', 'RSAmiss_PreResp_h', 'RSAhits_StimCue_h', 'RSAmiss_StimCue_h', 'RSAhits_StimResp_h', 'RSAmiss_StimResp_h');
    
elseif enoughH == 0% if not enough hippocampal SU
    save(['enoughH', subjID, '.mat'], 'enoughH');
end

%% cortex
cd ..\cortex
if enoughC == 1
    %% Hits (Cue/Cue) - 1
    [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_c, retSpikeNumber_cueLocked_c, 'hit', 0);
    RSAhits_CueCue_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c); % normalizes spikenumber and correlates enc+ret for RSA
    
    % Misses (Cue / Cue)
    if enoughMiss == 1 % only calculate this if there are enough misses
        [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_c, retSpikeNumber_cueLocked_c, 'miss', 0);
        RSAmiss_CueCue_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    else
        RSAmiss_CueCue_c = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Cue / Cue)';
    plotRSA(RSAhits_CueCue_c, RSAmiss_CueCue_c, comp, hitsIdx, missIdx);
    savefig(['RSA_cort_CueCue', subjID]);
    saveas(gcf,['RSA_cort_CueCue' subjID,'.png']);
    
    %% Hits (Cue/Resp) - 2
    [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_c, retSpikeNumber_respLocked_c, 'hit', 0);
    RSAhits_CueResp_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_c, retSpikeNumber_respLocked_c, 'miss', 0);
        RSAmiss_CueResp_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    else
        RSAmiss_CueResp_c = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Cue / Resp)';
    plotRSA(RSAhits_CueResp_c, RSAmiss_CueResp_c, comp, hitsIdx, missIdx);
    savefig(['RSA_cort_CueResp', subjID]);
    saveas(gcf,['RSA_cort_CueResp' subjID,'.png']);
    
    %% Hits (Resp/Cue) - 3
    [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_respLocked_c, retSpikeNumber_cueLocked_c, 'hit', 0);
    RSAhits_RespCue_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_respLocked_c, retSpikeNumber_cueLocked_c, 'miss', 0);
        RSAmiss_RespCue_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    else
        RSAmiss_RespCue_c = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Resp / Cue)';
    plotRSA(RSAhits_RespCue_c, RSAmiss_RespCue_c, comp, hitsIdx, missIdx);
    savefig(['RSA_cort_RespCue', subjID]);
    saveas(gcf,['RSA_cort_RespCue' subjID,'.png']);
    
    %% Hits (Resp / Resp) - 4
    [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_respLocked_c, retSpikeNumber_respLocked_c, 'hit', 0);
    RSAhits_RespResp_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_respLocked_c, retSpikeNumber_respLocked_c, 'miss', 0);
        RSAmiss_RespResp_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    else
        RSAmiss_RespResp_c = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Resp / Resp)';
    plotRSA(RSAhits_RespResp_c, RSAmiss_RespResp_c, comp, hitsIdx, missIdx);
    savefig(['RSA_cort_RespResp', subjID]);
    saveas(gcf,['RSA_cort_RespResp' subjID,'.png']);
    
    %% Hits (Pre / Cue) - 5
    [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_preTrial_c, retSpikeNumber_cueLocked_c, 'hit', 0);
    RSAhits_PreCue_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_preTrial_c, retSpikeNumber_cueLocked_c, 'miss', 0);
        RSAmiss_PreCue_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    else
        RSAmiss_PreCue_c = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Pre / Cue)';
    plotRSA(RSAhits_PreCue_c, RSAmiss_PreCue_c, comp, hitsIdx, missIdx);
    savefig(['RSA_cort_PreCue', subjID]);
    saveas(gcf,['RSA_cort_PreCue' subjID,'.png']);
    
    %% Hits (Pre / Resp) - 6
    [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_preTrial_c, retSpikeNumber_respLocked_c, 'hit', 0);
    RSAhits_PreResp_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    
    % PreResp Misses
    if enoughMiss == 1
        [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_preTrial_c, retSpikeNumber_respLocked_c, 'miss', 0);
        RSAmiss_PreResp_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    else
        RSAmiss_PreResp_c = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Pre / Resp)';
    plotRSA(RSAhits_PreResp_c, RSAmiss_PreResp_c, comp, hitsIdx, missIdx);
    savefig(['RSA_cort_PreCue', subjID]);
    saveas(gcf,['RSA_cort_PreCue' subjID,'.png']);
    
    %% Hits (Stim / Cue) - 7
    [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_c, retSpikeNumber_cueLocked_c, 'hit', 0);
    RSAhits_StimCue_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    
    % StimCue Misses
    if enoughMiss == 1
        [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_c, retSpikeNumber_cueLocked_c, 'miss', 0);
        RSAmiss_StimCue_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    else
        RSAmiss_StimCue_c = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Stim / Cue)';
    plotRSA(RSAhits_StimCue_c, RSAmiss_StimCue_c, comp, hitsIdx, missIdx);
    savefig(['RSA_cort_StimCue', subjID]);
    saveas(gcf,['RSA_cort_StimCue' subjID,'.png']);
    
    %% Hits (Stim / Resp) - 8
    [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_c, retSpikeNumber_respLocked_c, 'hit', 0);
    RSAhits_StimResp_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_c, tempCell_ret_c] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_c, retSpikeNumber_respLocked_c, 'miss', 0);
        RSAmiss_StimResp_c = normSpikeNumber(tempCell_enc_c, tempCell_ret_c);
    else
        RSAmiss_StimResp_c = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Stim / Resp)';
    plotRSA(RSAhits_StimResp_c, RSAmiss_StimResp_c, comp, hitsIdx, missIdx);
    savefig(['RSA_cort_StimResp', subjID]);
    saveas(gcf,['RSA_cort_StimResp' subjID,'.png']);
    
    %% save all RSAs
    save(['RSAmat_cort', subjID, '.mat'], 'RSAhits_CueCue_c', 'RSAmiss_CueCue_c', 'RSAhits_CueResp_c', 'RSAmiss_CueResp_c', 'RSAhits_RespCue_c', 'RSAmiss_RespCue_c', 'RSAhits_RespResp_c', 'RSAmiss_RespResp_c', 'RSAhits_PreCue_c', 'RSAmiss_PreCue_c', 'RSAhits_PreResp_c', 'RSAmiss_PreResp_c', 'RSAhits_StimCue_c', 'RSAmiss_StimCue_c', 'RSAhits_StimResp_c', 'RSAmiss_StimResp_c');
    
elseif enoughC == 0% if not enough cortical SU
    save(['enoughC', subjID, '.mat'], 'enoughC');
end

%% other
cd ..\other
if enoughO == 1
    %% Hits (Cue/Cue) - 1
    [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_o, retSpikeNumber_cueLocked_o, 'hit', 0);
    RSAhits_CueCue_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o); % normalizes spikenumber and correlates enc+ret for RSA
    
    % Misses (Cue / Cue)
    if enoughMiss == 1 % only calculate this if there are enough misses
        [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_o, retSpikeNumber_cueLocked_o, 'miss', 0);
        RSAmiss_CueCue_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    else
        RSAmiss_CueCue_o = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Cue / Cue)';
    plotRSA(RSAhits_CueCue_o, RSAmiss_CueCue_o, comp, hitsIdx, missIdx);
    savefig(['RSA_other_CueCue', subjID]);
    saveas(gcf,['RSA_other_CueCue' subjID,'.png']);
    
    %% Hits (Cue/Resp) - 2
    [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_o, retSpikeNumber_respLocked_o, 'hit', 0);
    RSAhits_CueResp_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_cueLocked_o, retSpikeNumber_respLocked_o, 'miss', 0);
        RSAmiss_CueResp_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    else
        RSAmiss_CueResp_o = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Cue / Resp)';
    plotRSA(RSAhits_CueResp_o, RSAmiss_CueResp_o, comp, hitsIdx, missIdx);
    savefig(['RSA_other_CueResp', subjID]);
    saveas(gcf,['RSA_other_CueResp' subjID,'.png']);
    
    %% Hits (Resp/Cue) - 3
    [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_respLocked_o, retSpikeNumber_cueLocked_o, 'hit', 0);
    RSAhits_RespCue_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_respLocked_o, retSpikeNumber_cueLocked_o, 'miss', 0);
        RSAmiss_RespCue_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    else
        RSAmiss_RespCue_o = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Resp / Cue)';
    plotRSA(RSAhits_RespCue_o, RSAmiss_RespCue_o, comp, hitsIdx, missIdx);
    savefig(['RSA_other_RespCue', subjID]);
    saveas(gcf,['RSA_other_RespCue' subjID,'.png']);
    
    %% Hits (Resp / Resp) - 4
    [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_respLocked_o, retSpikeNumber_respLocked_o, 'hit', 0);
    RSAhits_RespResp_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_respLocked_o, retSpikeNumber_respLocked_o, 'miss', 0);
        RSAmiss_RespResp_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    else
        RSAmiss_RespResp_o = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Resp / Resp)';
    plotRSA(RSAhits_RespResp_o, RSAmiss_RespResp_o, comp, hitsIdx, missIdx);
    savefig(['RSA_other_RespResp', subjID]);
    saveas(gcf,['RSA_other_RespResp' subjID,'.png']);
    
    %% Hits (Pre / Cue) - 5
    [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_preTrial_o, retSpikeNumber_cueLocked_o, 'hit', 0);
    RSAhits_PreCue_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_preTrial_o, retSpikeNumber_cueLocked_o, 'miss', 0);
        RSAmiss_PreCue_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    else
        RSAmiss_PreCue_o = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Pre / Cue)';
    plotRSA(RSAhits_PreCue_o, RSAmiss_PreCue_o, comp, hitsIdx, missIdx);
    savefig(['RSA_other_PreCue', subjID]);
    saveas(gcf,['RSA_other_PreCue' subjID,'.png']);
    
    %% Hits (Pre / Resp) - 6
    [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_preTrial_o, retSpikeNumber_respLocked_o, 'hit', 0);
    RSAhits_PreResp_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_preTrial_o, retSpikeNumber_respLocked_o, 'miss', 0);
        RSAmiss_PreResp_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    else
        RSAmiss_PreResp_o = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Pre / Resp)';
    plotRSA(RSAhits_PreResp_o, RSAmiss_PreResp_o, comp, hitsIdx, missIdx);
    savefig(['RSA_other_PreResp', subjID]);
    saveas(gcf,['RSA_other_PreResp' subjID,'.png']);
    
    %% Hits (Stim / Cue) - 7
    [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_o, retSpikeNumber_cueLocked_o, 'hit', 0);
    RSAhits_StimCue_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_o, retSpikeNumber_cueLocked_o, 'miss', 0);
        RSAmiss_StimCue_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    else
        RSAmiss_StimCue_o = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Stim / Cue)';
    plotRSA(RSAhits_StimCue_o, RSAmiss_StimCue_o, comp, hitsIdx, missIdx);
    savefig(['RSA_other_StimCue', subjID]);
    saveas(gcf,['RSA_other_StimCue' subjID,'.png']);
    
    %% Hits (Stim / Resp) - 8
    [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_o, retSpikeNumber_respLocked_o, 'hit', 0);
    RSAhits_StimResp_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    
    % respLocked Misses
    if enoughMiss == 1
        [tempCell_enc_o, tempCell_ret_o] = mk_tempCell(allTrials, encSpikeNumber_stimLocked_o, retSpikeNumber_respLocked_o, 'miss', 0);
        RSAmiss_StimResp_o = normSpikeNumber(tempCell_enc_o, tempCell_ret_o);
    else
        RSAmiss_StimResp_o = 0; % create a dummy variable to prevent error later when saving
    end
    
    % plotting
    comp = '(Stim / Resp)';
    plotRSA(RSAhits_StimResp_o, RSAmiss_StimResp_o, comp, hitsIdx, missIdx);
    savefig(['RSA_other_StimResp', subjID]);
    saveas(gcf,['RSA_other_StimResp' subjID,'.png']);
    
    %% save all RSAs
    save(['RSAmat_other', subjID, '.mat'], 'RSAhits_CueCue_o', 'RSAmiss_CueCue_o', 'RSAhits_CueResp_o', 'RSAmiss_CueResp_o', 'RSAhits_RespCue_o', 'RSAmiss_RespCue_o', 'RSAhits_RespResp_o', 'RSAmiss_RespResp_o', 'RSAhits_PreCue_o', 'RSAmiss_PreCue_o', 'RSAhits_PreResp_o', 'RSAmiss_PreResp_o', 'RSAhits_StimCue_o', 'RSAmiss_StimCue_o', 'RSAhits_StimResp_o', 'RSAmiss_StimResp_o');
    
elseif enoughO == 0% if not enough other SU
    save(['enoughO', subjID, '.mat'], 'enoughO');
end

%% RSA - Mask
% between category is 3
RSA_mask_hits = zeros(size(hitsIdx,1));
RSA_mask_hits = RSA_mask_hits+3;

% within category is 2
for x=1:size(hitsIdx,1)
    for y=1:size(hitsIdx,1)
        if x <=numHit(1) && y<=numHit(1)
            RSA_mask_hits(y,x)=2;
        elseif x>numHit(1) && x<=numHit(2) && y>numHit(1) && y<=numHit(2)
            RSA_mask_hits(y,x)=2;
        elseif x>numHit(2) && y>numHit(2)
            RSA_mask_hits(y,x)=2;
        end
    end
end

% main diagonal is 1
for x=1:size(hitsIdx,1)
    for y=1:size(hitsIdx,1)
        if x==y
            RSA_mask_hits(y,x)=1;
        end
    end
end

% Misses - Mask
% between category is 3
% i can probably just reuse the old RSA_mask variable
RSA_mask_miss=zeros(size(missIdx,1));
RSA_mask_miss=RSA_mask_miss+3;

% within category is 2
for x=1:size(missIdx,1)
    for y=1:size(missIdx,1)
        if x <=numMiss(1) && y<=numMiss(1)
            RSA_mask_miss(y,x)=2;
        elseif x>numMiss(1) && x<=numMiss(2) && y>numMiss(1) && y<=numMiss(2)
            RSA_mask_miss(y,x)=2;
        elseif x>numMiss(2) && y>numMiss(2)
            RSA_mask_miss(y,x)=2;
        end
    end
end

% main diagonal is 1
for x=1:size(missIdx,1)
    for y=1:size(missIdx,1)
        if x==y
            RSA_mask_miss(y,x)=1;
        end
    end
end

%% Cortex - Hits
if enoughC == 1
    
    if size(RSAhits_CueCue_c) == size(RSA_mask_hits)
        RSAmean_c(1,1) = mean(RSAhits_CueCue_c (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_c(1,2) = mean(RSAhits_CueCue_c (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_c(1,3) = mean(RSAhits_CueCue_c (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_CueResp_c) == size(RSA_mask_hits)
        RSAmean_c(2,1) = mean(RSAhits_CueResp_c (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_c(2,2) = mean(RSAhits_CueResp_c (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_c(2,3) = mean(RSAhits_CueResp_c (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_PreCue_c) == size(RSA_mask_hits)
        RSAmean_c(3,1) = mean(RSAhits_PreCue_c (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_c(3,2) = mean(RSAhits_PreCue_c (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_c(3,3) = mean(RSAhits_PreCue_c (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_PreResp_c) == size(RSA_mask_hits)
        RSAmean_c(4,1) = mean(RSAhits_PreResp_c (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_c(4,2) = mean(RSAhits_PreResp_c (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_c(4,3) = mean(RSAhits_PreResp_c (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_RespCue_c) == size(RSA_mask_hits)
        RSAmean_c(5,1) = mean(RSAhits_RespCue_c (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_c(5,2) = mean(RSAhits_RespCue_c (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_c(5,3) = mean(RSAhits_RespCue_c (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_RespResp_c) == size(RSA_mask_hits)
        RSAmean_c(6,1) = mean(RSAhits_RespResp_c (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_c(6,2) = mean(RSAhits_RespResp_c (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_c(6,3) = mean(RSAhits_RespResp_c (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_StimCue_c) == size(RSA_mask_hits)
        RSAmean_c(7,1) = mean(RSAhits_StimCue_c (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_c(7,2) = mean(RSAhits_StimCue_c (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_c(7,3) = mean(RSAhits_StimCue_c (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_StimResp_c) == size(RSA_mask_hits)
        RSAmean_c(8,1) = mean(RSAhits_StimResp_c (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_c(8,2) = mean(RSAhits_StimResp_c (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_c(8,3) = mean(RSAhits_StimResp_c (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    %% RSAmean_c - Miss
    if enoughMiss == 1
        
        if size(RSAmiss_CueCue_c) == size(RSA_mask_miss)
            RSAmean_c(1,4) = mean(RSAmiss_CueCue_c (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_c(1,5) = mean(RSAmiss_CueCue_c (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_c(1,6) = mean(RSAmiss_CueCue_c (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_CueResp_c) == size(RSA_mask_miss)
            RSAmean_c(2,4) = mean(RSAmiss_CueResp_c (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_c(2,5) = mean(RSAmiss_CueResp_c (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_c(2,6) = mean(RSAmiss_CueResp_c (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_PreCue_c) == size(RSA_mask_miss)
            RSAmean_c(3,4) = mean(RSAmiss_PreCue_c (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_c(3,5) = mean(RSAmiss_PreCue_c (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_c(3,6) = mean(RSAmiss_PreCue_c (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_PreResp_c) == size(RSA_mask_miss)
            RSAmean_c(4,4) = mean(RSAmiss_PreResp_c (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_c(4,5) = mean(RSAmiss_PreResp_c (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_c(4,6) = mean(RSAmiss_PreResp_c (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_RespCue_c) == size(RSA_mask_miss)
            RSAmean_c(5,4) = mean(RSAmiss_RespCue_c (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_c(5,5) = mean(RSAmiss_RespCue_c (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_c(5,6) = mean(RSAmiss_RespCue_c (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_RespResp_c) == size(RSA_mask_miss)
            RSAmean_c(6,4) = mean(RSAmiss_RespResp_c (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_c(6,5) = mean(RSAmiss_RespResp_c (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_c(6,6) = mean(RSAmiss_RespResp_c (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_StimCue_c) == size(RSA_mask_miss)
            RSAmean_c(7,4) = mean(RSAmiss_StimCue_c (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_c(7,5) = mean(RSAmiss_StimCue_c (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_c(7,6) = mean(RSAmiss_StimCue_c (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_StimResp_c) == size(RSA_mask_miss)
            RSAmean_c(8,4) = mean(RSAmiss_StimResp_c (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_c(8,5) = mean(RSAmiss_StimResp_c (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_c(8,6) = mean(RSAmiss_StimResp_c (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
    end
    
else
    RSAmean_c = [];
end

%% Hippocampus - Hits
if enoughH == 1
    
    if size(RSAhits_CueCue_h) == size(RSA_mask_hits)
        RSAmean_h(1,1) = mean(RSAhits_CueCue_h (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_h(1,2) = mean(RSAhits_CueCue_h (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_h(1,3) = mean(RSAhits_CueCue_h (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_CueResp_h) == size(RSA_mask_hits)
        RSAmean_h(2,1) = mean(RSAhits_CueResp_h (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_h(2,2) = mean(RSAhits_CueResp_h (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_h(2,3) = mean(RSAhits_CueResp_h (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_PreCue_h) == size(RSA_mask_hits)
        RSAmean_h(3,1) = mean(RSAhits_PreCue_h (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_h(3,2) = mean(RSAhits_PreCue_h (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_h(3,3) = mean(RSAhits_PreCue_h (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_PreResp_h) == size(RSA_mask_hits)
        RSAmean_h(4,1) = mean(RSAhits_PreResp_h (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_h(4,2) = mean(RSAhits_PreResp_h (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_h(4,3) = mean(RSAhits_PreResp_h (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_RespCue_h) == size(RSA_mask_hits)
        RSAmean_h(5,1) = mean(RSAhits_RespCue_h (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_h(5,2) = mean(RSAhits_RespCue_h (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_h(5,3) = mean(RSAhits_RespCue_h (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_RespResp_h) == size(RSA_mask_hits)
        RSAmean_h(6,1) = mean(RSAhits_RespResp_h (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_h(6,2) = mean(RSAhits_RespResp_h (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_h(6,3) = mean(RSAhits_RespResp_h (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_StimCue_h) == size(RSA_mask_hits)
        RSAmean_h(7,1) = mean(RSAhits_StimCue_h (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_h(7,2) = mean(RSAhits_StimCue_h (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_h(7,3) = mean(RSAhits_StimCue_h (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_StimResp_h) == size(RSA_mask_hits)
        RSAmean_h(8,1) = mean(RSAhits_StimResp_h (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_h(8,2) = mean(RSAhits_StimResp_h (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_h(8,3) = mean(RSAhits_StimResp_h (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    %% Hippocampus - Miss
    if enoughMiss == 1
        if size(RSAmiss_CueCue_h) == size(RSA_mask_miss)
            RSAmean_h(1,4) = mean(RSAmiss_CueCue_h (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_h(1,5) = mean(RSAmiss_CueCue_h (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_h(1,6) = mean(RSAmiss_CueCue_h (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_CueResp_h) == size(RSA_mask_miss)
            RSAmean_h(2,4) = mean(RSAmiss_CueResp_h (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_h(2,5) = mean(RSAmiss_CueResp_h (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_h(2,6) = mean(RSAmiss_CueResp_h (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_PreCue_h) == size(RSA_mask_miss)
            RSAmean_h(3,4) = mean(RSAmiss_PreCue_h (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_h(3,5) = mean(RSAmiss_PreCue_h (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_h(3,6) = mean(RSAmiss_PreCue_h (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_PreResp_h) == size(RSA_mask_miss)
            RSAmean_h(4,4) = mean(RSAmiss_PreResp_h (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_h(4,5) = mean(RSAmiss_PreResp_h (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_h(4,6) = mean(RSAmiss_PreResp_h (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_RespCue_h) == size(RSA_mask_miss)
            RSAmean_h(5,4) = mean(RSAmiss_RespCue_h (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_h(5,5) = mean(RSAmiss_RespCue_h (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_h(5,6) = mean(RSAmiss_RespCue_h (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_RespResp_h) == size(RSA_mask_miss)
            RSAmean_h(6,4) = mean(RSAmiss_RespResp_h (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_h(6,5) = mean(RSAmiss_RespResp_h (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_h(6,6) = mean(RSAmiss_RespResp_h (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_StimCue_h) == size(RSA_mask_miss)
            RSAmean_h(7,4) = mean(RSAmiss_StimCue_h (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_h(7,5) = mean(RSAmiss_StimCue_h (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_h(7,6) = mean(RSAmiss_StimCue_h (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_StimResp_h) == size(RSA_mask_miss)
            RSAmean_h(8,4) = mean(RSAmiss_StimResp_h (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_h(8,5) = mean(RSAmiss_StimResp_h (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_h(8,6) = mean(RSAmiss_StimResp_h (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
    end
    
else
    RSAmean_h = [];
end

%% Other - Hits
if enoughO == 1
    
    if size(RSAhits_CueCue_o) == size(RSA_mask_hits)
        RSAmean_o(1,1) = mean(RSAhits_CueCue_o (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_o(1,2) = mean(RSAhits_CueCue_o (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_o(1,3) = mean(RSAhits_CueCue_o (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_CueResp_o) == size(RSA_mask_hits)
        RSAmean_o(2,1) = mean(RSAhits_CueResp_o (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_o(2,2) = mean(RSAhits_CueResp_o (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_o(2,3) = mean(RSAhits_CueResp_o (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_PreCue_o) == size(RSA_mask_hits)
        RSAmean_o(3,1) = mean(RSAhits_PreCue_o (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_o(3,2) = mean(RSAhits_PreCue_o (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_o(3,3) = mean(RSAhits_PreCue_o (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_PreResp_o) == size(RSA_mask_hits)
        RSAmean_o(4,1) = mean(RSAhits_PreResp_o (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_o(4,2) = mean(RSAhits_PreResp_o (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_o(4,3) = mean(RSAhits_PreResp_o (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_RespCue_o) == size(RSA_mask_hits)
        RSAmean_o(5,1) = mean(RSAhits_RespCue_o (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_o(5,2) = mean(RSAhits_RespCue_o (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_o(5,3) = mean(RSAhits_RespCue_o (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_RespResp_o) == size(RSA_mask_hits)
        RSAmean_o(6,1) = mean(RSAhits_RespResp_o (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_o(6,2) = mean(RSAhits_RespResp_o (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_o(6,3) = mean(RSAhits_RespResp_o (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_StimCue_o) == size(RSA_mask_hits)
        RSAmean_o(7,1) = mean(RSAhits_StimCue_o (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_o(7,2) = mean(RSAhits_StimCue_o (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_o(7,3) = mean(RSAhits_StimCue_o (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    if size(RSAhits_StimResp_o) == size(RSA_mask_hits)
        RSAmean_o(8,1) = mean(RSAhits_StimResp_o (RSA_mask_hits == 1) ); % withinTrial
        RSAmean_o(8,2) = mean(RSAhits_StimResp_o (RSA_mask_hits == 2) ); % withinCateg
        RSAmean_o(8,3) = mean(RSAhits_StimResp_o (RSA_mask_hits == 3) ); % betweenCate
    else
        error('no match between RSA matrix and RSA mask');
    end
    
    %% Other - Miss
    if enoughMiss == 1
        if size(RSAmiss_CueCue_o) == size(RSA_mask_miss)
            RSAmean_o(1,4) = mean(RSAmiss_CueCue_o (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_o(1,5) = mean(RSAmiss_CueCue_o (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_o(1,6) = mean(RSAmiss_CueCue_o (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_CueResp_o) == size(RSA_mask_miss)
            RSAmean_o(2,4) = mean(RSAmiss_CueResp_o (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_o(2,5) = mean(RSAmiss_CueResp_o (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_o(2,6) = mean(RSAmiss_CueResp_o (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_PreCue_o) == size(RSA_mask_miss)
            RSAmean_o(3,4) = mean(RSAmiss_PreCue_o (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_o(3,5) = mean(RSAmiss_PreCue_o (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_o(3,6) = mean(RSAmiss_PreCue_o (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_PreResp_o) == size(RSA_mask_miss)
            RSAmean_o(4,4) = mean(RSAmiss_PreResp_o (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_o(4,5) = mean(RSAmiss_PreResp_o (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_o(4,6) = mean(RSAmiss_PreResp_o (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_RespCue_o) == size(RSA_mask_miss)
            RSAmean_o(5,4) = mean(RSAmiss_RespCue_o (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_o(5,5) = mean(RSAmiss_RespCue_o (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_o(5,6) = mean(RSAmiss_RespCue_o (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_RespResp_o) == size(RSA_mask_miss)
            RSAmean_o(6,4) = mean(RSAmiss_RespResp_o (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_o(6,5) = mean(RSAmiss_RespResp_o (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_o(6,6) = mean(RSAmiss_RespResp_o (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_StimCue_o) == size(RSA_mask_miss)
            RSAmean_o(7,4) = mean(RSAmiss_StimCue_o (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_o(7,5) = mean(RSAmiss_StimCue_o (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_o(7,6) = mean(RSAmiss_StimCue_o (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
        
        if size(RSAmiss_StimResp_o) == size(RSA_mask_miss)
            RSAmean_o(8,4) = mean(RSAmiss_StimResp_o (RSA_mask_miss == 1) ); % withinTrial
            RSAmean_o(8,5) = mean(RSAmiss_StimResp_o (RSA_mask_miss == 2) ); % withinCateg
            RSAmean_o(8,6) = mean(RSAmiss_StimResp_o (RSA_mask_miss == 3) ); % betweenCate
        else
            error('no match between RSA matrix and RSA mask');
        end
    end
    
else
    RSAmean_o = [];
end

% this output is better for SPSS
if enoughMiss == 0
    suffix = [1 1; 2 1; 3 1];
else
    suffix = [1,1;2,1;3,1;1,0;2,0;3,0];
end

if ~isempty(RSAmean_c)
    RSAmean_c = RSAmean_c';
    RSAmean_c = [RSAmean_c suffix];
end

if ~isempty(RSAmean_h)
    RSAmean_h = RSAmean_h';
    RSAmean_h = [RSAmean_h suffix];
end

if ~isempty(RSAmean_o)
    RSAmean_o = RSAmean_o';
    RSAmean_o = [RSAmean_o suffix];
end

cd ../..
save(['RSAmean_', subjID,'.mat'], 'RSAmean_c', 'RSAmean_h', 'RSAmean_o');

end

clear allSubj
for sesh = 1:40
allSubj{sesh,1} = strcat('P99_S', num2str(sesh));
end