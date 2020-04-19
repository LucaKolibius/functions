%% firing onset between index neurons during encoding and retrieval
cd X:\Luca\indexSUvisu\encRet_p95_th1
clear; clc
load('globalTest.mat')
encVec_allSU(encVec_allSU>1 | encVec_allSU<0) = [];
retVec_allSU(retVec_allSU>1 | retVec_allSU<0) = [];

% gaussian kernel
mlength = [-0.075:0.0002:0.075];
mSigma  = 0.02; % 20ms
mSigma  = 0.01;
mKernel = normpdf(mlength,0,mSigma);
mKernel = mKernel/max(mKernel); % normalize peak to 1

% encoding
dt        = [0:0.001:1]; % from 0s to +1s in steps of 1ms
encHist   = histogram(encVec_allSU,dt);
encHist   = encHist.BinCounts;
encConv   = conv(mKernel, encHist);
% get rid of edges
encConv(1:375)       = [];
encConv(end-374:end) = [];

% retrieval
retHist = histogram(retVec_allSU,dt);
retHist = retHist.BinCounts;
retConv = conv(mKernel,retHist);
% get rid of edges
retConv(1:375)       = [];
retConv(end-374:end) = [];

% plot convolved series
figure; subplot(211); hold on;
plot(encConv, 'linew', 3); plot(retConv,'linew',3);
ylabel('Smoothed Occurences')
mAx = gca;
mAx.YAxis.FontWeight = 'bold';
mAx.XAxis.FontWeight = 'bold';
mAx.FontSize = 16;
box on
legend({'Encoding', 'Retrieval'}, 'FontSize',16, 'FontWeight', 'bold');
legend boxoff
xlabel('Time after cue onset (ms)')

subplot(212); hold on;
[a,lags] = xcorr(encConv,retConv, 'coeff');
plot(lags,a, 'linew', 3)
ylim([0 1])
plot([0 0], [get(gca, 'YLim')], '--', 'Color', [0.6 0.6 0.6], 'linew', 3)
ylabel('Correlation');
xlabel('Lag in ms')
mAx = gca;
mAx.YAxis.FontWeight = 'bold';
mAx.XAxis.FontWeight = 'bold';
mAx.FontSize = 16;

figure
[lagsA, A] = xcorr(encConv, 'coeff');
subplot(211); hold on;
plot(A, lagsA)
subplot(212)
[lagsA, A] = xcorr(retConv, 'coeff');
plot(A, lagsA)


%% permutation testing
nperm = 10000;
xcorPerm = zeros(nperm, size(encConv,2)*2-1); % preallocation (10000x1999 zero double)
for perm = 1:nperm
    % possibly I dont have to save every single permutation but just
    % calculate the permutated bins + log the corrcoeff of that to then
    % use prctile-95 to find a signiifcance threshold
    encHist_perm = encHist(randperm(size(encHist,2)));  % maybe dont even save this
    retHist_perm = retHist(randperm(size(retHist,2)));
    
    encConv_perm   = conv(mKernel, encHist_perm);
    % get rid of edges
    encConv_perm(1:375)       = [];
    encConv_perm(end-374:end) = [];
    
    retConv_perm = conv(mKernel,retHist_perm);
    % get rid of edges
    retConv_perm(1:375)       = [];
    retConv_perm(end-374:end) = [];
    
    [xcorP, ~] = xcorr(encConv_perm, retConv_perm, 'coeff');
    xcorPerm   = [xcorPerm; xcorP];
    
%     subplot(311)
%     plot(encConv_perm)
%     subplot(312)
%     plot(retConv_perm)
%     subplot(313)
%     plot(xcorP)
%     pause(3)
end

permTH = prctile(xcorPerm, 95);
figure; plot(permTH);
%% sin2 comes after sin1
figure; subplot(211)
t = [-0.5:0.01:0.5];
f = 3;
sin1 = sin(2*pi*f*t);
plot(t,sin1, 'linew',3)

sin2 = sin(2*pi*f*t-0.5);
hold on; plot(t,sin2, 'linew',3);
legend({'Sin1', 'Sin2'}, 'FontSize',16, 'FontWeight','bold');

subplot(212)
[a,lag] = xcorr(sin1,sin2, 'coeff');
plot(lag,a, 'linew',3)

% for xcorr(c,d) a peak that is below 0 means that "c" comes before "d"