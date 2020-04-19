%% permutation testing
% all SU time series between encoding and retrieval
clear
close
fetch_allSubj
permCor = zeros(377, 10000);
hitMiss = 'hit';
encRet = 'enc'; % 'enc' or 'ret'
plvl = 0.01; % 5%
numperm = 10000; % number of permutations
maxLag = 1; % maximum lag for the autocorrelation; in trials
rownum = 0; % counter for the total number of SU
empSignAR = 0; % counter how many auto correlations are signficant
for sbj = 1 : size(allSubj,1)
    subjID = allSubj{sbj};
    disp(subjID);
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
    cd advancedAnalysis\elecLoc\spikenumbers
    
    % load tables for enc and ret
    abc = dir('spikeNumber_hipp_*');
    load(abc.name, 'encSpikeNumber_cueLocked_h', 'retSpikeNumber_cueLocked_h')
    
    allTrials = size(encSpikeNumber_cueLocked_h,2);
    if size(encSpikeNumber_cueLocked_h,1)>2
        [tempCell_enc_h, tempCell_ret_h] = mk_tempCell2(allTrials, encSpikeNumber_cueLocked_h, retSpikeNumber_cueLocked_h, hitMiss); % extracts the absolute spikenumbers for encoding and retrieval. CAREFUL: only hits and ordered by ff/pp/fp
    else
        continue
    end
    
    if size(tempCell_enc_h,2)>1 % in case there are no misses
        [~ , normENC, normRET] = normSpikeNumber(tempCell_enc_h, tempCell_ret_h); % normalizes spikenumber and correlates enc+ret for RSA
    else
        continue
    end
    subjID(regexp(subjID,'_')) = '-';
    
    if strcmp(encRet,'enc')
        orgDat = normENC;
    elseif strcmp(encRet, 'ret')
        orgDat = normRET;
    end
    
    cd X:\Luca\data\SUvisu
    % loop over all single units in that session
    for su = 1 : size(orgDat,1)
        x = orgDat(su,:);
        [orgAR, lag] = xcorr(x,x, maxLag, 'coeff');
        orgAR(maxLag+1) = [];
        orgRT = sum(orgAR);
                rownum = rownum+1;

        % permutation
        for prm = 1:numperm
            permTS = orgDat(1, randperm(size(orgDat,2)));
            [r, ~] = xcorr(permTS,permTS, maxLag, 'coeff');
            r(maxLag+1) = []; % get rid of lag-0;
            rT = sum(r);
            permCor(rownum,prm) = rT;
        end
        
        % we need a threshold now
        % first get rid of 1/-1 rows and NaNs
        if sum(isnan(permCor(rownum,:))) == numperm || (sum(permCor(rownum,:)==1 | permCor(rownum,:)==-1) == numperm)% if all permutations of that SU are nans or either -1 or 1
            permCor(rownum,:) = [];
            rownum = rownum-1;
            disp(sprintf('problem at %.0f', rownum));
        else % if I keep the permutation for that SU, calculate the threshold
            forTH = sort(permCor(rownum,:));
            thresh(rownum) = forTH(numperm-(plvl*numperm));
        end
                
        if orgRT>=thresh(rownum)
            empSignAR = empSignAR+1;
        end
    end
end

%% extract randomly one permutation from each SU
numsign = zeros(1,10000);
for ih = 1:numperm
    permdxASU = randperm(numperm,size(permCor,1))'; % creates 377 (which is the number of SU) random numbers from 1 to 10000 (number of permutations)
    temp = permCor([1:size(permCor,1)], permdxASU);
    permvec(:,ih) = temp(logical(eye(size(permCor,1)))); % a vector of randomly drawn values from all SU (one value per SU)
    numsign(ih) = sum(permvec(:,ih)>=thresh'); % the number of randomly drawn h0 correlations that are signficiantly above threshold
end

%%
%% retrieval reordered to presentation sequence
[~, ~, ~, ~, ~, ~, ~, orda] = loadLogs(p2d);
normRET_orda = normRET(:,orda);
signNiv(su) = 1-(sum(orgTS>shufR)/numperm);


%% playing around
clear temp
temp = sort(numsign,'ascend')
temp(9500)
empSignAR
temp(end)
