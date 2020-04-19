% using this script I deduced that a plvl of .01 should be used for testing
% whether the encoding and retrieval TS correlate significantly. Likewise a
% plvl of .01 should be used for determining a suprathreshold dp (dot
% product) indicating the trial to which the single unit is sensitive to
close
figure(26)
clear
hitsAll     = [];
falseposAll = [];
signTScorAll = [];

for rep = 1:100
for ix = 1:11 % steps of signal/noise
pause(2)

% random numbers
% encTS = randi(20, [1,40]);
% retTS = randi(20, [1,40]);
encTS = 1+randn([1 40]);
retTS = 1+randn([1 40]);

% signal
signal = [1:0.2:3];

% add signal to trial 10 encoding and retrieval
encTS(10) = encTS(10)+signal(ix);
retTS(10) = retTS(10)+signal(ix)/3*2;

% same with trial 12
encTS(12) = encTS(10)+signal(ix);
retTS(12) = retTS(10)+signal(ix)/2;

% draw TS for enc and ret
subplot(311)
hold off
plot(encTS); hold on; plot(retTS);

% calc and draw D-P
dp = encTS .* retTS;
subplot(312)
hold off
plot(dp, 'g');
hold on
title('Dot-Product of Trial Series Above')

% calculate correlation between the original trial series (encoding + retrieval)
temp = corrcoef(encTS, retTS);
origCor(ix) = temp(2);

% permutate to get distribution
for nperm = 1:10000
permDex = randperm(40,40); % index for reshuffling retTS
permRet(:,nperm) = retTS(permDex); % reshuffled ret TS

permDex = randperm(40,40); % repeat for encoding
permEnc(:,nperm) = encTS(permDex);

temp = corrcoef(encTS, permRet(:,nperm)); % calculate new correlation
permCor(nperm) = temp(2); % save new correlation

newDP(nperm,:) = encTS .* permRet(:,nperm)';
end
thresh(ix)   = prctile(permCor,99); % cor = .99
threshDP(ix) = prctile(newDP(:),99); % dp = .99

plot([1 40], [threshDP(ix) threshDP(ix)], 'g--');

subplot(311)
hold on
% ret in red, enc in blue
threshRet(ix) = prctile(permRet(:), 99); % threshold for significant retrieval trials
threshEnc(ix) = prctile(permEnc(:), 99); % threshold for significant encoding  trials
plot([1 40], [threshRet(ix) threshRet(ix)], 'r--');
plot([1 40], [threshEnc(ix) threshEnc(ix)], 'b--');


numsignDP(ix) = sum(dp>=threshDP(ix));
hits(ix) = double(dp(10)>=threshDP(ix)) + double(dp(12)>=threshDP(ix));
falsepos(ix) = numsignDP(ix)-hits(ix);

title(sprintf('Number signficant: %.0f', numsignDP(ix)));
subplot(313)
xlim([1 11])
plot(thresh,'r');
hold on
plot(origCor, 'b');
end

signTScor      = origCor > thresh;
signTScorAll   = [signTScorAll; signTScor];
hitsAll        = [hitsAll; hits];
falseposAll    = [falseposAll; falsepos];
end
% save('vars99cor99dp.mat', 'signTScorAll', 'hitsAll', 'falseposAll');

%% visualize results
clear
clc
% load('vars99cor99dp.mat');

% all TS considered hitrate
allTS_hitrate = mean(hitsAll);
hold on
allTS_fp = mean(falseposAll);

% get out significant TS
signTS_hits = [];
signTS_fp   = [];
for tp = 1:11
  signTS_hits = [signTS_hits,  mean(hitsAll(logical(signTScorAll(:,tp)),tp))];
  signTS_fp   = [signTS_fp, mean(falseposAll(logical(signTScorAll(:,tp)),tp))];
end

figure
hold on
plot(allTS_hitrate, 'r', 'linew',5);
plot(allTS_fp, 'g', 'linew',5);
plot(signTS_hits, 'r', 'linew', 1);
plot(signTS_fp,'g', 'linew',1);

ylim([0 2.5]);