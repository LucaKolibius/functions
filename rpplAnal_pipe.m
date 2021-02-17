%% NEW VERSION FOR MICRO LFP WITHOUT SPKINT THAT ONLY CONSIDERES BUNDLES IN WHICH I HAVE HIPPOCAMPAL UNITS
function rpplAnal_pipe

% calculate the ripples on each bundle (take the wire with the highest ripple activity)
% then look over the whole bundle which trials are indexed (rppl_stats2 l.35)
% only use the first SU in each bundle (rppl_stats2 l.25)

whereAmI(0)
global prePath;
addpath([prePath, 'Luca\functions']);
addpath([prePath, 'Luca\toolboxes\fieldtrip-20200310']); ft_defaults;
addpath([prePath, 'Luca\toolboxes\fieldtrip-20200603']); ft_defaults;
rmpath(genpath('\\analyse4.psy.gla.ac.uk\project0310\RDS\Common\mcode'));
% lfpDir = dir([prePath, 'Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_noSPKINT.mat']);  % CHANGED TO NO SPK INT
lfpDir = dir([prePath, 'Luca\data\microLFP\sub-*_onlyMicroLFP_RAW_1000DS_*.mat']);  % spkInt in laptop and noSpkInt on local
load([prePath, 'Luca\data\allSbj\allSpksHZ.mat'], 'allSpks')

allSUPow.idx  = [];
allSUPow.ndx  = [];
allSUPow.dff  = [];
allFreqRes    = [];

for spk = 1 : length(allSpks)
    
        if any(isnan(allSpks(spk).idxTrlSingHi)) % this basically saves some computation time. would be evenbetter if i would only take the first su of each bundle
            continue
        end
    
    %% GET: bidsID + sesh
    bidsID  = allSpks(spk).bidsID;
    sesh    = allSpks(spk).sesh;
    %     allBund = {allSpks.bundlename};
    curBund = allSpks(spk).bundlename;
    %     curWire = [curBund, num2str(allSpks(spk).favChan)];
    %     favChan = allSpks(spk).favChan(45:50);
    idxTrl  = allSpks(spk).idxTrlSingHi;     % WHICH TRIALS DOES THAT WIRE INDEX?
    encTrig = round(allSpks(spk).encTrigger(allSpks(spk).hitsIdx,[1 3])*1000);
    
    %% LOAD IN THE LFP-DATA
    % There is no indendent LFP for sub-1007_S1b (it is in the same
    % recording as S1, it's just another session)
    load([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_noSPKINT.mat'], 'data');
    
    % LAPTOP AND NO SPK INT
    %         abc = dir([lfpDir(1).folder, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_onlyMicroLFP_RAW_1000DS_*.mat']);
    %         load([abc.folder, filesep, abc.name], 'data');
    
    %% REDEFINE TRIALS ACCORDING TO encTrig
    %     cfg          = [ ];
    %     cfg.trl      = [encTrig(:,1)-1000-100-200 encTrig(:,2)+100 zeros(size(encTrig,1))];
    %     microLFPlong = ft_redefinetrial(cfg, data);
    %     microLFPlong = rmfield(microLFPlong, 'trialinfo');
    
    cfg          = [ ];
    cfg.trl      = [encTrig(:,1) encTrig(:,2) zeros(size(encTrig,1))];
    microLFP     = ft_redefinetrial(cfg, data);
    microLFP     = rmfield(microLFP, 'trialinfo');
    
    %% DEMEAN & ORTHOGONALIZE
    cfg        = [];
    cfg.demean = 'yes';
    microLFP   = ft_preprocessing(cfg, microLFP);
    microLFP   = orthogonVec(microLFP);           % orthogonalize
    
    %% SELECT CHANNEL THAT REFLECTS SU INPUT (see spkPhase_encRet)
    %     cfg          = [];
    %     cfg.channel  = curWire;
    %     microLFP     = ft_selectdata(cfg, microLFP); % select
    
    %% CHANNEL SELECTION COULD BE A DIFFERENT WIRE FOR EACH FREQUENCY
    allBund = cellfun(@(x) x(1:end-1), microLFP.label, 'un', 0);
    bundIdx = strcmp(allBund, curBund);
    
    cfg         = [];
    cfg.channel = microLFP.label(bundIdx);
    microLFP    = ft_selectdata(cfg, microLFP);
    
    %% DETECT RIPPLES
    [ripple, bndLab, peaks, staEnd, favChan] = calcRppl (microLFP);
%     allSpks(spk).rpplNum = ripple.number;
%     allSpks(spk).rpplLen = ripple.length;
%     allSpks(spk).rpplDen = ripple.density;
    
    
% %%    GET RAW AVERAGE RIPPLE SHAPE
%     rppl.filt = [];
%     rppl.lfp  = [];
%     for trl = 1: size(peaks,1)
%         %% CONTINUE IF TRIAL HAS NO RIPPLES
%         if isempty(peaks{trl,1})
%             continue
%         end
%         
%         rpNum = size(peaks{trl,1},2);
%         for rip = 1 : rpNum
%             
%             
%             LFP    = microLFP.trial{trl}(favChan,:);
%             filtLFP = ft_preproc_bandpassfilter(microLFP.trial{trl}(favChan,:), 1000, [80 140], 3*fix(1000/80)+1, 'fir', 'twopass');
%             
% %             mStart = peaks{trl,1}(rip)-50; if mStart<1; mStart = 1; end
% %             mEnd   = peaks{trl,2}(rip)+50; if mEnd>size(microLFP.trial{trl},2); mEnd = size(microLFP.trial{trl},2); end
%             mStart = peaks{trl}(rip)-500; if mStart<1; continue; end
%             mEnd   = peaks{trl}(rip)+500; if mEnd>size(microLFP.trial{trl},2); continue; end
%             
%             LFPrip     = LFP(mStart:mEnd);
%             filtLFPrip = filtLFP(mStart:mEnd);
%             
% %             % zero pad
% %             zeropad = (400-size(LFPrip,2))/2;
% %             switch round(zeropad) == zeropad % is zeropad a full number
% %                 
% %                 case 0 % has decimals
% %                     LFPrip     = [ zeros(1,floor(zeropad)+1) LFPrip        zeros(1,floor(zeropad)) ];
% %                     filtLFPrip = [ zeros(1,floor(zeropad)+1) filtLFPrip    zeros(1,floor(zeropad)) ];
% %                 case 1  % no decimals
% %                     LFPrip     = [ zeros(1,floor(zeropad))   LFPrip        zeros(1,floor(zeropad)) ];
% %                     filtLFPrip = [ zeros(1,floor(zeropad))   filtLFPrip    zeros(1,floor(zeropad)) ];
% %             
% %             end
% %           
% %             if size(LFPrip,2) ~= 400
% %                 error('lfp snippet is not 300 long')
% %             end
%             
%             
%             rppl(spk).lfp  = [ rppl(spk).lfp;  LFPrip  ];
%             rppl(spk).filt = [ rppl(spk).filt; filtLFPrip ];
%         end
%         
%     end

%% get spike density estimate over trial
rpplTms = zeros(length(staEnd),20*1000); 
for trl = 1 : length(staEnd)
    if isempty(staEnd{trl,1})
        continue
    end
    
    mstart = staEnd{trl,1};
    mend   = staEnd{trl,2};
    
    for rip = 1 : size(mstart,2)
        rpplTms(trl,mstart(rip):mend(rip)) = 1;
    end  
    
end
allSpks(spk).rpplTms = rpplTms;

    
    
% end % END OF TRIAL LOOP

end % END OF SU LOOP

%%
%     % CALCULATE TRIAL RIPPLE POWER (80-140hz)
%     [idxTrlPow, ndxTrlPow, goodTrl, freqRes] = calcRipplPow (microLFP, idxTrl, encTrig, favChan);
%
%     if sum(idxTrl(logical(goodTrl))) > 0
%
%         allSUPow.idx  = [ allSUPow.idx  {idxTrlPow}  ];
%         allSUPow.ndx  = [ allSUPow.ndx  {ndxTrlPow}  ];
%
%     end

%     allFreqRes = [allFreqRes; freqRes];


% save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_rppls80to140.mat', 'allSpks', '-v7.3');
% save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\avgRppl.mat', 'rppl', '-v7.3'); % ==> rpplVisu_average_raw_filt
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\rpplDens_time.mat', 'allSpks', '-v7.3');

end % END OF FUNCTION