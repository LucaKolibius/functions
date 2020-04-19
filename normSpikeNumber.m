% normalize spikenumber & correlate enc + ret
function [RSAmat, encMat, retMat]= normSpikeNumber(tempCell_enc, tempCell_ret)
clusterMean = [];
clusterStd = [];
for i=1:size(tempCell_enc,1) % for number of SU
    clusterMean(i,1) = mean([tempCell_enc{i,:} tempCell_ret{i,:}]);
    clusterStd(i,1) = std([tempCell_enc{i,:} tempCell_ret{i,:}]);
end


for ia = 1:size(tempCell_enc,2)
    for ib = 1:size(tempCell_enc,1)
        encMat(ib,ia) = (tempCell_enc{ib,ia} - clusterMean(ib,1)) / clusterStd(ib,1);
        retMat(ib,ia) = (tempCell_ret{ib,ia} - clusterMean(ib,1)) / clusterStd(ib,1);
    end
end

% get rid of clusters that do not fire in any trial
delClust = clusterMean ~= 0; % indexing
encMat = encMat(delClust,:); % delete respective cluster
retMat = retMat(delClust,:); % delete respective cluster

if size(encMat,1) < 2 % if it is only one or no SU left
    RSAmat = nan( [size(tempCell_enc,2), size(tempCell_enc,2)] ); % size depends on the number of miss trials (mXm NaN matrix where m is the number of misses)
else % otherwise
    RSAmat=[];
    for iEnc = 1 : size(encMat,2)
        for iRet = 1 : size(retMat,2)
            mCorrelation = corrcoef(encMat(:,iEnc),retMat(:,iRet));
            RSAmat(iRet,iEnc) = mCorrelation(2);
        end
    end
    
end

%% delete the ones that are 0 and keep the nans?
% that will just lead to a nan matrix...
% wouldve worked as well. this way I create a NaN matrix if there are not
% enough SU anymore
