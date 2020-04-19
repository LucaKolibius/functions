% Subsequent Memory Effect (more firing during encoding of later remembered trials)
% Extract spikes within -1s to 2s cue locked for retrieval and -1 to +5 for encoding

function SME(subjID)

mSubject = subjID(1:end-3);
mSession = subjID(end-1:end);

try
    cd Z:/Luca/data
catch
    cd /media/ldk898/rds-share/Luca/data
end

if regexp(mSubject, '_')
    mSubject(end) = [];
    mSession = ['S', mSession];
end

cd(mSubject)
cd(mSession)
abc = dir;
cd(abc(3).name)
cd advancedAnalysis
mkdir('SME');
cd SME
SMEdir = cd;
cd ../postXCrej/
loadVar = dir('tableTimestamps_postXC_*');
load(loadVar.name);
cd ../manualRej/
loadVar = dir('clNames*');
load(loadVar.name);

% load logfile
cd ../..
p2d = cd;
p2d(end+1) = '\';
[~, hitsIdx, missIdx, ~, ~, ~, encTrigger] = loadLogs(p2d);

% gaussian kernel
mlength = [-0.075:0.0002:0.075];
mSigma  = 0.02; % 20ms
mKernel = normpdf(mlength,0,mSigma);
mKernel = mKernel/max(mKernel); % normalize peak to 1

% Encoding
clCounter   = 0;
gausConvENC = [];
dt          = [-1:0.0002:5]; % steps of 0.2ms
dt_bl       = [-1:0.0002:0];
for ia=1:size(tableTimestamps_postXC,2)
    mVarname = tableTimestamps_postXC.Properties.VariableNames{1,ia};
    if ismember(mVarname, clNames)
        clCounter     = clCounter + 1;
        clusterSpikes = table2array(tableTimestamps_postXC{1,ia})'; % load in all spiketimes from cluster mVarname
        clusterSpikes = clusterSpikes/32000; % samples to sec
        % segment trials of this SU into trials
        for ib = 1:size(encTrigger,1) % number of trials
            mBinary = clusterSpikes(clusterSpikes >= (encTrigger(ib, 1) -1) ...
                & clusterSpikes <= (encTrigger(ib,1) + 5))  - encTrigger(ib,1);
            mBinary_bl = clusterSpikes(clusterSpikes >= (encTrigger(ib,1)-1) ...
                & clusterSpikes <= (encTrigger(ib,1) + 0)) - encTrigger(ib,1);
            
            % do hist and convolution for raw
            mhist  = hist(mBinary,dt);
            myTS   = conv(mKernel, mhist);
            
            % do hist and convolution for bl
            mhist_bl = hist(mBinary_bl,dt_bl);
            myTS_bl  = conv(mKernel, mhist_bl);
            
            % correct for bins
            myTS([1:750, end-750:end])    = [];
            myTS_bl([1:750, end-750:end]) = [];
            
            % substract mean(bl) from raw
            subsBL     = mean(myTS_bl);
            myTS_blcor = myTS - subsBL;

            % write in output
            gausConvENC(clCounter, :, ib) = myTS_blcor; % SU x time series x trials
        end
    end
end

gausConvENC_hits = gausConvENC(:, : ,hitsIdx);
gausConvENC_miss = gausConvENC(:, :, missIdx);

gausConvENC_hits = sum(gausConvENC_hits,1);
gausConvENC_miss = sum(gausConvENC_miss,1);

gausConvENC_hits = mean(gausConvENC_hits,3);
gausConvENC_miss = mean(gausConvENC_miss,3);

% % e voila
% gausConvENC_hits = mean(mean(gausConvENC_hits,1),3);
% gausConvENC_miss = mean(mean(gausConvENC_miss,1),3);

% saving variable james ryan
cd(SMEdir)
save(['SME_', subjID, '.mat'], 'gausConvENC_hits', 'gausConvENC_miss');


% reshaped_hits = mean(reshape(gausConv_hits, [], 29250), 1);
% test = reshape(gausConv_hits, [], 29250);
% 
% test2= reshape(gausConv_hits, 29250, []);
% 
% mean(gausConv_hits_m)
% mean(reshaped_hits)
% mean(mean(test2))
% 
% 
% reshaped_miss = mean(reshape(gausConv_miss, [], 29250), 1);


  
% %% #4 average across trials (one for all, one for hits, one for misses)
% % #5 baseline correction
% % encoding - hits
% temp_bl = gausConvENC_bl; % data used for baseline
% temp_bl = mean(mean(temp_bl));
% rawTS = gausConvENC(hitsIdx,:); % uncorrected time series
% 
% %blCorrected = rawTS - temp_bl; % correct rawTS with BL
% gausConvENC_hit = mean(blCorrected) - temp_bl;
% 
% % encoding - miss
% rawTS = gausConvENC(missIdx,:);
% blCorrected = rawTS - temp_bl;
% gausConvENC_miss = mean(blCorrected);
% 
% % same for retrieval
% % retrieval - hits
% temp_bl = gausConvRET_bl;
% temp_bl = mean(mean(temp_bl));
% 
% rawTS = gausConvRET(hitsIdx,:);
% blCorrected = rawTS - temp_bl;
% gausConvRET_hit = mean(blCorrected);
% 
% % retrieval - miss
% rawTS = gausConvRET(missIdx,:);
% blCorrected = rawTS - temp_bl;
% gausConvRET_miss = mean(blCorrected);
end
