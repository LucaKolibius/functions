% normalize spikenumber
% works with double input instead of cell
% discards SU that do not fire in any trial
function [encMat, retMat]= normSpikeNumber_double(tempCell_enc, tempCell_ret)
clusterMean = [];
clusterMeanEnc = [];
clusterMeanRet = [];
clusterStd = [];
clusterStdEnc = [];
clusterStdRet = [];

for i=1:size(tempCell_enc,1) % for number of SU
    clusterMean(i,1) = mean([tempCell_enc(i,:) tempCell_ret(i,:)]);
    clusterMeanEnc(i,1) = mean(tempCell_enc(i,:));
    clusterMeanRet(i,1) = mean(tempCell_ret(i,:));
    
    clusterStd(i,1) = std([tempCell_enc(i,:) tempCell_ret(i,:)]);
    clusterStdEnc(i,1) = std(tempCell_enc(i,:));
    clusterStdRet(i,1) = std(tempCell_ret(i,:));
end

% for ia = 1:size(tempCell_enc,2)
%     for ib = 1:size(tempCell_enc,1)
%         encMat(ib,ia) = (tempCell_enc(ib,ia) - clusterMean(ib,1)) / clusterStd(ib,1);
%         retMat(ib,ia) = (tempCell_ret(ib,ia) - clusterMean(ib,1)) / clusterStd(ib,1);
%     end
% end

    for ia = 1:size(tempCell_enc,1)
        encMat(ia,:) = (tempCell_enc(ia,:) - clusterMeanEnc(ia)) / clusterStdEnc(ia);
        retMat(ia,:) = (tempCell_ret(ia,:) - clusterMeanRet(ia)) / clusterStdRet(ia);
    end

% get rid of clusters that do not fire in any trial
delClust = clusterMean ~= 0; % indexing
encMat = encMat(delClust,:); % delete respective cluster
retMat = retMat(delClust,:); % delete respective cluster