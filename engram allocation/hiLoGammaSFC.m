clear

%% SERVER
% load('/analyse/Project0309/Luca/data/allSbj/allPhsPowDat_34hz_70hz.mat', 'allPhsDat');
% addpath('/analyse/Project0309/Luca/toolboxes/CircStat2012a')

%% PC
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allPhsPowDat_34hz_70hz.mat', 'allPhsDat');

%% START SCRIPT
foi = 34:2:70;
fromTo = 0:30:360;
fromTo = deg2rad(fromTo);
count = 1;

for su = 1 : size(allPhsDat,2)
    
    lineLength = fprintf('%d SU out of %d SU done (%.2f%%).\n', su, length(allPhsDat), su/length(allPhsDat)*100);

    if isempty(allPhsDat(su).encIdx)
        continue
    end
    
    % AT LEAST 20 SPIKES REQUIRED
    numDat                 =  [allPhsDat(su).encIdx(1,:) allPhsDat(su).retIdx(1,:)];
    numDat(isnan(numDat))  =  [];
    numDat                 =  size(numDat,2);
   
    if numDat < 20
        continue
    end
    
    phsDatSU = [allPhsDat(su).encIdx allPhsDat(su).retIdx];

    %% HIGHEST LOW FREQ SFC
    zval = [];
    for freq = 2:6
    
        dat             = phsDatSU(freq,:);
        dat(isnan(dat)) = [];
        [~,zval(freq)]  = circ_rtest(dat);
    
    end
    [~, posLow] = max(zval);

    %% HIGHEST HIGH FREQ SFC
    zval = [];
    for freq = 14:19
        
        dat             = phsDatSU(freq,:);
        dat(isnan(dat)) = [];
         [~,zval(freq)]  = circ_rtest(dat);
        
    end
    [~, posHigh] = max(zval);

    phsLw = phsDatSU(posLow,:)  + pi;
    phsHi = phsDatSU(posHigh,:) + pi;
    
    %% LOW TO HIGH
    % FIND LOW FREQUENCY PEAK
    polHand        = polarhistogram(phsLw, fromTo, 'Normalization', 'probability');
    [~, maxPos]    = max(polHand.Values);
    BOI            = [maxPos maxPos+1]; % BIN OF INTEREST
     
    % ORIGINAL HIGH POLAR PLOT
    polHand        = polarhistogram(phsHi, fromTo, 'Normalization', 'probability');
    hiValOrg       = polHand.Values;
    
    % TRIMMED HIGH POLAR PLOT
    idx            = phsHi >= fromTo(BOI(1)) & phsHi <= fromTo(BOI(2));
    phsHiTrim      = phsHi;
    phsHiTrim(idx) = [];
    
    polHand        = polarhistogram(phsHiTrim, fromTo, 'normalization', 'probability');
        hiValTrim      = polHand.Values;
    
    diffValHi_emp      = hiValTrim-hiValOrg;
    
    
    % PERMUTATION
    nperm = 10000;
    diffValHi_perm = zeros(nperm, size(diffValHi_emp,2));
    for perm = 1:nperm
        idxPerm                = idx(randperm(size(idx,2)));
        phsHiTrimPerm          = phsHi;
        phsHiTrimPerm(idxPerm) = [];
        polHand                = polarhistogram(phsHiTrimPerm, fromTo, 'normalization', 'probability');
        hiValTrimPerm          = polHand.Values;
        
        diffValHi_perm(perm,:) = hiValTrimPerm-hiValOrg;
        
    end

    % NORMALIZE THROUGH PERMUTATION
    diffValHi_zscr = diffValHi_emp-mean(diffValHi_perm,1);
    diffValHi_zscr = diffValHi_zscr ./ std(diffValHi_perm, 0, 1);
    
    low2high(count,:) = diffValHi_zscr;
    
        %% REPEAT FOR HIGH LOW
    % FIND LOW FREQUENCY PEAK
    polHand        = polarhistogram(phsHi, fromTo, 'Normalization', 'probability');
    [~, maxPos]    = max(polHand.Values);
    BOI            = [maxPos maxPos+1]; % BIN OF INTEREST
     
    % ORIGINAL HIGH POLAR PLOT
    polHand        = polarhistogram(phsLw, fromTo, 'Normalization', 'probability');
    lwValOrg       = polHand.Values;
    
    % TRIMMED HIGH POLAR PLOT
    idx            = phsLw >= fromTo(BOI(1)) & phsLw <= fromTo(BOI(2));
    phsLwTrim      = phsLw;
    phsLwTrim(idx) = [];
    
    polHand        = polarhistogram(phsLwTrim, fromTo, 'normalization', 'probability');
    lwValTrim      = polHand.Values;
    
    diffValLw_emp  = lwValTrim-lwValOrg;
    
    
    % PERMUTATION
    nperm = 10000;
        diffValLw_perm = zeros(nperm, size(diffValLw_emp,2));
    for perm = 1:nperm
        idxPerm                = idx(randperm(size(idx,2))); % randomly take other phase values
        phsLwTrimPerm          = phsLw;
        phsLwTrimPerm(idxPerm) = []; % trim
        polHand                = polarhistogram(phsLwTrimPerm, fromTo, 'normalization', 'probability');
        lwValTrimPerm          = polHand.Values;
        
        diffValLw_perm(perm,:) = lwValTrimPerm-lwValOrg;
        
    end

    % NORMALIZE THROUGH PERMUTATION
    diffValLw_zscr = diffValLw_emp -   mean(diffValLw_perm,1);
    diffValLw_zscr = diffValLw_zscr ./ std(diffValLw_perm, 0, 1);
    
    high2low(count,:) = diffValLw_zscr;
    
   count = count + 1; 
   
   fprintf(repmat('\b',1,lineLength-0));

end