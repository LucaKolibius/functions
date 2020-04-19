% nsess=40;
%
% for n=1:nsess
%     matsiz=round(rand(1,1)*10+40);
%     mat(n).temp=rand(matsiz,matsiz)+eye(matsiz,matsiz);
%     matres(:,:,n)=imresize(mat(n).temp,[1000 1000]);
% end
%
% AVG=mean(matres,3);
%
% figure;imagesc(AVG);

function Average4RSA(allSubj)
%% make SQmask
% preallocation
% fetch_allSubj
for it = 1 : size(allSubj,1)
    subjID = allSubj{it};
    
    [numHit, numMiss] = get_FfFpPp_hits_miss(subjID); % this is the absolute number of ff/fp/pp instead of adding the previous number of ff/fp/pp to the next one!
    
    %     % alternative writing
    %     for i = 1:9
    %     Sq_new(i).hi = i*ones(numHit (ceil(i/3)), numHit (mod(i+2,3)+1));
    %     end
    
    % HITS - mask
    % actual size
    Sq1_hits = 1*ones(numHit(1),numHit(1));
    Sq2_hits = 2*ones(numHit(1),numHit(2));
    Sq3_hits = 3*ones(numHit(1),numHit(3));
    Sq4_hits = 4*ones(numHit(2),numHit(1));
    Sq5_hits = 5*ones(numHit(2),numHit(2));
    Sq6_hits = 6*ones(numHit(2),numHit(3));
    Sq7_hits = 7*ones(numHit(3),numHit(1));
    Sq8_hits = 8*ones(numHit(3),numHit(2));
    Sq9_hits = 9*ones(numHit(3),numHit(3));
    
    % put the individual rectangles back to a mask RSA
    SQ(it).mask_hits         =  [Sq1_hits Sq2_hits Sq3_hits; Sq4_hits Sq5_hits Sq6_hits; Sq7_hits Sq8_hits Sq9_hits];
    SQ(it).size_hits(:,:)    =  [size(Sq1_hits)' size(Sq2_hits)' size(Sq3_hits)' size(Sq4_hits)' size(Sq5_hits)' size(Sq6_hits)' size(Sq7_hits)' size(Sq8_hits)' size(Sq9_hits)']; % size of each rectangle
    
    % MISS
    Sq1_miss = 1*ones(numMiss(1),numMiss(1));
    Sq2_miss = 2*ones(numMiss(1),numMiss(2));
    Sq3_miss = 3*ones(numMiss(1),numMiss(3));
    Sq4_miss = 4*ones(numMiss(2),numMiss(1));
    Sq5_miss = 5*ones(numMiss(2),numMiss(2));
    Sq6_miss = 6*ones(numMiss(2),numMiss(3));
    Sq7_miss = 7*ones(numMiss(3),numMiss(1));
    Sq8_miss = 8*ones(numMiss(3),numMiss(2));
    Sq9_miss = 9*ones(numMiss(3),numMiss(3));
    
    SQ(it).mask_miss         =  [Sq1_miss Sq2_miss Sq3_miss; Sq4_miss Sq5_miss Sq6_miss; Sq7_miss Sq8_miss Sq9_miss];
    SQ(it).size_miss(:,:)    =  [size(Sq1_miss)' size(Sq2_miss)' size(Sq3_miss)' size(Sq4_miss)' size(Sq5_miss)' size(Sq6_miss)' size(Sq7_miss)' size(Sq8_miss)' size(Sq9_miss)'];
    
end
cd ../../
save('SQ.mat', 'SQ'); % this variable will have the mask RSA and the size of each rectangle for hits and misses

% calculating the maximum dimensions of each square (I know its actually a
% rectangle, fuck off) over all sessions. (:,:,1) are hits and (:,:,2) are
% misses
allVal = vertcat(SQ.size_hits);
val1 = allVal(1:2:size(allVal,1),:);
val2 = allVal(2:2:size(allVal,1),:);
max1 = max(val1,[],1);
max2 = max(val2,[],1);
allMax(:,:,1) = [max1; max2];

allVal = vertcat(SQ.size_miss);
val1 = allVal(1:2:size(allVal,1),:);
val2 = allVal(2:2:size(allVal,1),:);
max1 = max(val1,[],1);
max2 = max(val2,[],1);
allMax(:,:,2) = [max1; max2];

% manually change it to 32x32 plus a mold
allMax(:,:,1) = [38 38 38 32 32 32 38 38 38; 38 32 38 38 32 38 38 32 38];

clearvars -except SQ allSubj allMax

%% apply that to all 9x2 RSAs
% cycle through all subjects and sessions
for sbj = 1 : size(allSubj,1)
    % cd to correct folder
    try
        cd X:/Luca/data
    catch
        cd /media/ldk898/rds-share/Luca/data
    end
    
    subjID = allSubj{sbj};
    mSubject = subjID(1:end-3);
    mSession = subjID(end-1:end);
    
    % if the session name is called 1b ("one-b") then this line prevents an error during cd
    mSubject(regexp(mSubject,'_')) = [];
    if isempty(regexp(mSession,'S', 'ONCE'))
        mSession = ['S', mSession];
    end
    
    cd(mSubject)
    cd(mSession)
    abc = dir; cd(abc(3).name)
    cd advancedAnalysis\RSA\oldSize
    
    %% cycle through ROI
    dROI = cd;
    for roi = 1:3
        if roi == 1
            cd other
            if ~exist('rzRSAhits_other', 'var')
                rzRSAhits_other = struct('CueCue', {[]}, 'CueResp',  {[]},'PreCue', {[]},'PreResp', {[]},'RespCue', {[]},'RespResp', {[]},'StimCue',{[]}, 'StimResp',{[]});
                rzRSAmiss_other = struct('CueCue', {[]}, 'CueResp',  {[]},'PreCue', {[]},'PreResp', {[]},'RespCue', {[]},'RespResp', {[]},'StimCue',{[]}, 'StimResp',{[]});
            end
        elseif roi == 2
            cd cortex
            if ~exist('rzRSAhits_cort', 'var')
                rzRSAhits_cort = struct('CueCue', {[]}, 'CueResp',  {[]},'PreCue', {[]},'PreResp', {[]},'RespCue', {[]},'RespResp', {[]},'StimCue',{[]}, 'StimResp',{[]});
                rzRSAmiss_cort = struct('CueCue', {[]}, 'CueResp',  {[]},'PreCue', {[]},'PreResp', {[]},'RespCue', {[]},'RespResp', {[]},'StimCue',{[]}, 'StimResp',{[]});
            end
        elseif roi == 3
            cd hippocampus
            if ~exist('rzRSAhits_hipp', 'var')
                rzRSAhits_hipp = struct('CueCue', {[]}, 'CueResp',  {[]},'PreCue', {[]},'PreResp', {[]},'RespCue', {[]},'RespResp', {[]},'StimCue',{[]}, 'StimResp',{[]});
                rzRSAmiss_hipp = struct('CueCue', {[]}, 'CueResp',  {[]},'PreCue', {[]},'PreResp', {[]},'RespCue', {[]},'RespResp', {[]},'StimCue',{[]}, 'StimResp',{[]});
            end
        end
        
        % load all RSAs
        clear RSA*
        
        enoughX = dir('enough*.mat');
        if isempty(enoughX)
            abc = dir('RSAmat*'); load(abc.name)
        else
            cd(dROI)
            continue
        end
        % cycle through RSAs
        allRSA = who('RSA*');
        
        for ix = 1 : size(allRSA,1)
            oldRSA = [];
            oldRSA = eval(allRSA{ix});
            
            % HITS
            if ~isempty(regexp(allRSA{ix},'hits','ONCE'))
                for sqx = 1:9 % 9 squares
                    oldSize = SQ(sbj).size_hits(:,sqx); % get the old size of square number sqx
                    
                    
                    newSize = allMax(:,sqx,1); % whole column / sq1 / hits
                    
                    % extract oldSq from RSA
                    mySquare = oldRSA(SQ(sbj).mask_hits==sqx); % 1 for sq1
                    % oldSq(:,:,sqx) = reshape(temp, [actualSize(:,1)]');
                    resSQ(sqx).oldSq = reshape(mySquare, oldSize');
                    
                    % resize
                    % newSq(:,:,sqx) = imresize(oldSq, [goalSize(:,1)]');
                    resSQ(sqx).newSq = imresize(resSQ(sqx).oldSq, newSize');
                end
            elseif ~isempty(regexp(allRSA{ix},'miss','ONCE'))
                % MISS
                for sqx = 1:9 % 9 squares
                    % sq1
                    oldSize = SQ(sbj).size_miss(:,sqx);
                    newSize = allMax(:,sqx,2); % whole column / sqx / miss
                    
                    % extract oldSq from RSA
                    mySquare = oldRSA(SQ(sbj).mask_miss==sqx); % 1 for sq1
                    resSQ(sqx).oldSq = reshape(mySquare, oldSize'); % reshape into rectangle
                    
                    if isempty(mySquare) % if there is no square put NaNs into it
                        resSQ(sqx).newSq = NaN(newSize');
                        continue
                    else
                        
                        % resize
                        resSQ(sqx).newSq = imresize(resSQ(sqx).oldSq, newSize');
                    end
                end
            end
            resRSA = [resSQ(1).newSq resSQ(2).newSq resSQ(3).newSq; resSQ(4).newSq resSQ(5).newSq resSQ(6).newSq; resSQ(7).newSq resSQ(8).newSq resSQ(9).newSq];
            
            %% Sort resRSA into the correct variable
            temp = allRSA{ix};
            if ~isempty(regexp(temp,'_h','ONCE')) % HIPPOCAMPUS
                if ~isempty(regexp(temp,'hits','ONCE')) % HITS
                    if ~isempty(regexp(temp,'CueCue','ONCE'))
                        rzRSAhits_hipp(sbj).CueCue = resRSA;
                    elseif ~isempty(regexp(temp,'CueResp','ONCE'))
                        rzRSAhits_hipp(sbj).CueResp = resRSA;
                    elseif ~isempty(regexp(temp,'PreCue','ONCE'))
                        rzRSAhits_hipp(sbj).PreCue = resRSA;
                    elseif ~isempty(regexp(temp,'PreResp','ONCE'))
                        rzRSAhits_hipp(sbj).PreResp = resRSA;
                    elseif ~isempty(regexp(temp,'RespCue','ONCE'))
                        rzRSAhits_hipp(sbj).RespCue = resRSA;
                    elseif ~isempty(regexp(temp,'RespResp','ONCE'))
                        rzRSAhits_hipp(sbj).RespResp = resRSA;
                    elseif ~isempty(regexp(temp,'StimCue','ONCE'))
                        rzRSAhits_hipp(sbj).StimCue = resRSA;
                    elseif ~isempty(regexp(temp,'StimResp','ONCE'))
                        rzRSAhits_hipp(sbj).StimResp = resRSA;
                    end
                elseif ~isempty(regexp(temp,'miss','ONCE')) %% MISS
                    if ~isempty(regexp(temp,'CueCue','ONCE'))
                        rzRSAmiss_hipp(sbj).CueCue = resRSA;
                    elseif ~isempty(regexp(temp,'CueResp','ONCE'))
                        rzRSAmiss_hipp(sbj).CueResp = resRSA;
                    elseif ~isempty(regexp(temp,'PreCue','ONCE'))
                        rzRSAmiss_hipp(sbj).PreCue = resRSA;
                    elseif ~isempty(regexp(temp,'PreResp','ONCE'))
                        rzRSAmiss_hipp(sbj).PreResp = resRSA;
                    elseif ~isempty(regexp(temp,'RespCue','ONCE'))
                        rzRSAmiss_hipp(sbj).RespCue = resRSA;
                    elseif ~isempty(regexp(temp,'RespResp','ONCE'))
                        rzRSAmiss_hipp(sbj).RespResp = resRSA;
                    elseif ~isempty(regexp(temp,'StimCue','ONCE'))
                        rzRSAmiss_hipp(sbj).StimCue = resRSA;
                    elseif ~isempty(regexp(temp,'StimResp','ONCE'))
                        rzRSAmiss_hipp(sbj).StimResp = resRSA;
                    end
                end
            elseif ~isempty(regexp(temp,'_c','ONCE')) %% CORTEX
                if ~isempty(regexp(temp,'hits','ONCE')) % HITS
                    if ~isempty(regexp(temp,'CueCue','ONCE'))
                        rzRSAhits_cort(sbj).CueCue = resRSA;
                    elseif ~isempty(regexp(temp,'CueResp','ONCE'))
                        rzRSAhits_cort(sbj).CueResp = resRSA;
                    elseif ~isempty(regexp(temp,'PreCue','ONCE'))
                        rzRSAhits_cort(sbj).PreCue = resRSA;
                    elseif ~isempty(regexp(temp,'PreResp','ONCE'))
                        rzRSAhits_cort(sbj).PreResp = resRSA;
                    elseif ~isempty(regexp(temp,'RespCue','ONCE'))
                        rzRSAhits_cort(sbj).RespCue = resRSA;
                    elseif ~isempty(regexp(temp,'RespResp','ONCE'))
                        rzRSAhits_cort(sbj).RespResp = resRSA;
                    elseif ~isempty(regexp(temp,'StimCue','ONCE'))
                        rzRSAhits_cort(sbj).StimCue = resRSA;
                    elseif ~isempty(regexp(temp,'StimResp','ONCE'))
                        rzRSAhits_cort(sbj).StimResp = resRSA;
                    end
                elseif ~isempty(regexp(temp,'miss','ONCE')) %% MISS
                    if ~isempty(regexp(temp,'CueCue','ONCE'))
                        rzRSAmiss_cort(sbj).CueCue = resRSA;
                    elseif ~isempty(regexp(temp,'CueResp','ONCE'))
                        rzRSAmiss_cort(sbj).CueResp = resRSA;
                    elseif ~isempty(regexp(temp,'PreCue','ONCE'))
                        rzRSAmiss_cort(sbj).PreCue = resRSA;
                    elseif ~isempty(regexp(temp,'PreResp','ONCE'))
                        rzRSAmiss_cort(sbj).PreResp = resRSA;
                    elseif ~isempty(regexp(temp,'RespCue','ONCE'))
                        rzRSAmiss_cort(sbj).RespCue = resRSA;
                    elseif ~isempty(regexp(temp,'RespResp','ONCE'))
                        rzRSAmiss_cort(sbj).RespResp = resRSA;
                    elseif ~isempty(regexp(temp,'StimCue','ONCE'))
                        rzRSAmiss_cort(sbj).StimCue = resRSA;
                    elseif ~isempty(regexp(temp,'StimResp','ONCE'))
                        rzRSAmiss_cort(sbj).StimResp = resRSA;
                    end
                end
            elseif ~isempty(regexp(temp,'_o','ONCE')) %% OTHER
                if ~isempty(regexp(temp,'hits','ONCE')) % HITS
                    if ~isempty(regexp(temp,'CueCue','ONCE'))
                        rzRSAhits_other(sbj).CueCue = resRSA;
                    elseif ~isempty(regexp(temp,'CueResp','ONCE'))
                        rzRSAhits_other(sbj).CueResp = resRSA;
                    elseif ~isempty(regexp(temp,'PreCue','ONCE'))
                        rzRSAhits_other(sbj).PreCue = resRSA;
                    elseif ~isempty(regexp(temp,'PreResp','ONCE'))
                        rzRSAhits_other(sbj).PreResp = resRSA;
                    elseif ~isempty(regexp(temp,'RespCue','ONCE'))
                        rzRSAhits_other(sbj).RespCue = resRSA;
                    elseif ~isempty(regexp(temp,'RespResp','ONCE'))
                        rzRSAhits_other(sbj).RespResp = resRSA;
                    elseif ~isempty(regexp(temp,'StimCue','ONCE'))
                        rzRSAhits_other(sbj).StimCue = resRSA;
                    elseif ~isempty(regexp(temp,'StimResp','ONCE'))
                        rzRSAhits_other(sbj).StimResp = resRSA;
                    end
                elseif ~isempty(regexp(temp,'miss','ONCE')) %% MISS
                    if ~isempty(regexp(temp,'CueCue','ONCE'))
                        rzRSAmiss_other(sbj).CueCue = resRSA;
                    elseif ~isempty(regexp(temp,'CueResp','ONCE'))
                        rzRSAmiss_other(sbj).CueResp = resRSA;
                    elseif ~isempty(regexp(temp,'PreCue','ONCE'))
                        rzRSAmiss_other(sbj).PreCue = resRSA;
                    elseif ~isempty(regexp(temp,'PreResp','ONCE'))
                        rzRSAmiss_other(sbj).PreResp = resRSA;
                    elseif ~isempty(regexp(temp,'RespCue','ONCE'))
                        rzRSAmiss_other(sbj).RespCue = resRSA;
                    elseif ~isempty(regexp(temp,'RespResp','ONCE'))
                        rzRSAmiss_other(sbj).RespResp = resRSA;
                    elseif ~isempty(regexp(temp,'StimCue','ONCE'))
                        rzRSAmiss_other(sbj).StimCue = resRSA;
                    elseif ~isempty(regexp(temp,'StimResp','ONCE'))
                        rzRSAmiss_other(sbj).StimResp = resRSA;
                    end
                end
            end
            
            %% change directory to corresponding ROI in folder newSize
            mcd = cd;
            indx = regexp(mcd,'oldSize\')+8;
            cROI = mcd(indx:end);
            cd ../../newSize
            if ~exist('hippocampus','file'); mkdir('hippocampus'); end
            if ~exist('cortex','file'); mkdir('cortex'); end
            if ~exist('other','file'); mkdir('other'); end
            cd(cROI)
            
            % save resized RSA aka rzRSA (8 for hits and 8 for misses = 16 variables per Session)
            save(['rz', allRSA{ix}], 'resRSA');
            cd(mcd)
        end
        %% end of roi
        cd(dROI) % go back to RSA folder
        
    end
end

cd X:\Luca\data\analysis\RSA_visu\
save(['rzRSA_hipp',  subjID, '.mat'], 'rzRSAhits_hipp',  'rzRSAmiss_hipp');
save(['rzRSA_cort',  subjID, '.mat'], 'rzRSAhits_cort',  'rzRSAmiss_cort');
save(['rzRSA_other', subjID, '.mat'], 'rzRSAhits_other', 'rzRSAmiss_other');


%% Average the resized RSA
% HIPP HIT
sizeHits = size(rzRSAhits_hipp(1).CueCue);

mAverage = zeros([sizeHits,8]);
dividedby = 0;
for it = 1:size(rzRSAhits_hipp,2)
    if isempty(rzRSAhits_hipp(it).CueCue)
        continue
    end
    dividedby = dividedby + 1;
    mAverage(:,:,1) = mAverage(:,:,1) + rzRSAhits_hipp(it).CueCue;
    mAverage(:,:,2) = mAverage(:,:,2) + rzRSAhits_hipp(it).CueResp;
    mAverage(:,:,3) = mAverage(:,:,3) + rzRSAhits_hipp(it).PreCue;
    mAverage(:,:,4) = mAverage(:,:,4) + rzRSAhits_hipp(it).PreResp;
    mAverage(:,:,5) = mAverage(:,:,5) + rzRSAhits_hipp(it).RespCue;
    mAverage(:,:,6) = mAverage(:,:,6) + rzRSAhits_hipp(it).RespResp;
    mAverage(:,:,7) = mAverage(:,:,7) + rzRSAhits_hipp(it).StimCue;
    mAverage(:,:,8) = mAverage(:,:,8) + rzRSAhits_hipp(it).StimResp;
end

figH = figure(1);
hold on
subplot(4,2,1)
CueCue_h = imagesc(mAverage(:,:,1)/dividedby);
title('Cue / Cue')
colorbar

subplot(4,2,2)
CueResp_h = imagesc(mAverage(:,:,2)/dividedby);
title('Cue / Resp')
colorbar

subplot(4,2,3)
PreCue_h = imagesc(mAverage(:,:,3)/dividedby);
title('Pre / Cue')
colorbar

subplot(4,2,4)
PreResp_h = imagesc(mAverage(:,:,4)/dividedby);
title('Pre / Resp')
colorbar

subplot(4,2,5)
RespCue_h = imagesc(mAverage(:,:,5)/dividedby);
title('Resp / Cue')
colorbar

subplot(4,2,6)
RespResp_h = imagesc(mAverage(:,:,6)/dividedby);
title('Resp / Resp')
colorbar

subplot(4,2,7)
StimCue_h = imagesc(mAverage(:,:,7)/dividedby);
title('Stim / Cue')
colorbar

subplot(4,2,8)
StimResp_h = imagesc(mAverage(:,:,8)/dividedby);
title('Stim / Resp')
colorbar

%% CORT HIT
mAverage = zeros([sizeHits,8]);
dividedby = 0;
for it = 1:size(rzRSAhits_cort,2)
    if isempty(rzRSAhits_cort(it).CueCue)
        continue
    end
    dividedby = dividedby + 1;
    mAverage(:,:,1) = mAverage(:,:,1) + rzRSAhits_cort(it).CueCue;
    mAverage(:,:,2) = mAverage(:,:,2) + rzRSAhits_cort(it).CueResp;
    mAverage(:,:,3) = mAverage(:,:,3) + rzRSAhits_cort(it).PreCue;
    mAverage(:,:,4) = mAverage(:,:,4) + rzRSAhits_cort(it).PreResp;
    mAverage(:,:,5) = mAverage(:,:,5) + rzRSAhits_cort(it).RespCue;
    mAverage(:,:,6) = mAverage(:,:,6) + rzRSAhits_cort(it).RespResp;
    mAverage(:,:,7) = mAverage(:,:,7) + rzRSAhits_cort(it).StimCue;
    mAverage(:,:,8) = mAverage(:,:,8) + rzRSAhits_cort(it).StimResp;
end

figH = figure(1);
hold on
subplot(4,2,1)
CueCue_h = imagesc(mAverage(:,:,1)/dividedby);
title('Cue / Cue')
colorbar

subplot(4,2,2)
CueResp_h = imagesc(mAverage(:,:,2)/dividedby);
title('Cue / Resp')
colorbar

subplot(4,2,3)
PreCue_h = imagesc(mAverage(:,:,3)/dividedby);
title('Pre / Cue')
colorbar

subplot(4,2,4)
PreResp_h = imagesc(mAverage(:,:,4)/dividedby);
title('Pre / Resp')
colorbar

subplot(4,2,5)
RespCue_h = imagesc(mAverage(:,:,5)/dividedby);
title('Resp / Cue')

subplot(4,2,6)
RespResp_h = imagesc(mAverage(:,:,6)/dividedby);
title('Resp / Resp')
colorbar

subplot(4,2,7)
StimCue_h = imagesc(mAverage(:,:,7)/dividedby);
title('Stim / Cue')
colorbar

subplot(4,2,8)
StimResp_h = imagesc(mAverage(:,:,8)/dividedby);
title('Stim / Resp')
colorbar

%% HIPP MISS
sizeMiss = size(rzRSAmiss_hipp(1).CueCue);
mAverage = zeros([sizeMiss,8]);
dividedby = 0;
for it = 1:size(rzRSAmiss_hipp,2)
    if isempty(rzRSAmiss_hipp(it).CueCue)
        continue
    end
    dividedby = dividedby + 1;
    mAverage(:,:,1) = mAverage(:,:,1) + rzRSAmiss_hipp(it).CueCue;
    mAverage(:,:,2) = mAverage(:,:,2) + rzRSAmiss_hipp(it).CueResp;
    mAverage(:,:,3) = mAverage(:,:,3) + rzRSAmiss_hipp(it).PreCue;
    mAverage(:,:,4) = mAverage(:,:,4) + rzRSAmiss_hipp(it).PreResp;
    mAverage(:,:,5) = mAverage(:,:,5) + rzRSAmiss_hipp(it).RespCue;
    mAverage(:,:,6) = mAverage(:,:,6) + rzRSAmiss_hipp(it).RespResp;
    mAverage(:,:,7) = mAverage(:,:,7) + rzRSAmiss_hipp(it).StimCue;
    mAverage(:,:,8) = mAverage(:,:,8) + rzRSAmiss_hipp(it).StimResp;
end

figH = figure(1);
hold on
subplot(4,2,1)
CueCue_h = imagesc(mAverage(:,:,1)/dividedby);
title('Cue / Cue')

subplot(4,2,2)
CueResp_h = imagesc(mAverage(:,:,2)/dividedby);
title('Cue / Resp')

subplot(4,2,3)
PreCue_h = imagesc(mAverage(:,:,3)/dividedby);
title('Pre / Cue')

subplot(4,2,4)
PreResp_h = imagesc(mAverage(:,:,4)/dividedby);
title('Pre / Resp')

subplot(4,2,5)
RespCue_h = imagesc(mAverage(:,:,5)/dividedby);
title('Resp / Cue')

subplot(4,2,6)
RespResp_h = imagesc(mAverage(:,:,6)/dividedby);
title('Resp / Resp')

subplot(4,2,7)
StimCue_h = imagesc(mAverage(:,:,7)/dividedby);
title('Stim / Cue')

subplot(4,2,8)
StimResp_h = imagesc(mAverage(:,:,8)/dividedby);
title('Stim / Resp')

%% CORT MISS
mAverage = zeros([sizeMiss,8]);
dividedby = 0;
for it = 1:size(rzRSAmiss_cort,2)
    if isempty(rzRSAmiss_cort(it).CueCue)
        continue
    end
    dividedby = dividedby + 1;
    mAverage(:,:,1) = mAverage(:,:,1) + rzRSAmiss_cort(it).CueCue;
    mAverage(:,:,2) = mAverage(:,:,2) + rzRSAmiss_cort(it).CueResp;
    mAverage(:,:,3) = mAverage(:,:,3) + rzRSAmiss_cort(it).PreCue;
    mAverage(:,:,4) = mAverage(:,:,4) + rzRSAmiss_cort(it).PreResp;
    mAverage(:,:,5) = mAverage(:,:,5) + rzRSAmiss_cort(it).RespCue;
    mAverage(:,:,6) = mAverage(:,:,6) + rzRSAmiss_cort(it).RespResp;
    mAverage(:,:,7) = mAverage(:,:,7) + rzRSAmiss_cort(it).StimCue;
    mAverage(:,:,8) = mAverage(:,:,8) + rzRSAmiss_cort(it).StimResp;
end

figH = figure(1);
hold on
subplot(4,2,1)
CueCue_h = imagesc(mAverage(:,:,1)/dividedby);
title('Cue / Cue')

subplot(4,2,2)
CueResp_h = imagesc(mAverage(:,:,2)/dividedby);
title('Cue / Resp')

subplot(4,2,3)
PreCue_h = imagesc(mAverage(:,:,3)/dividedby);
title('Pre / Cue')

subplot(4,2,4)
PreResp_h = imagesc(mAverage(:,:,4)/dividedby);
title('Pre / Resp')

subplot(4,2,5)
RespCue_h = imagesc(mAverage(:,:,5)/dividedby);
title('Resp / Cue')

subplot(4,2,6)
RespResp_h = imagesc(mAverage(:,:,6)/dividedby);
title('Resp / Resp')

subplot(4,2,7)
StimCue_h = imagesc(mAverage(:,:,7)/dividedby);
title('Stim / Cue')

subplot(4,2,8)
StimResp_h = imagesc(mAverage(:,:,8)/dividedby);
title('Stim / Resp')


end % end of function

% respresp & cueresp
% do a RSA that is not reordered?