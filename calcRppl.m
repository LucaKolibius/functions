function [rpplRec, rpplsWire, bndLab] = calcRppl (data, bidsID, sesh, artFold)

chanLen   =  length( data.time{1} ) ./ 1000;  % timepoint length divided by sampling freq 1024 Hz gives channellength in seconds
rppl_exp  =  chanLen * (0.023 + 0.00545);           % number of ripples in white noise of same length

% reject artefact wires
load([artFold, filesep, bidsID, '_', sesh, '_micro_th-8.mat'], 'del_sampleinfo', 'delIx')

% NAMES OF THE BUNDLES
bndLab = unique(cellfun(@(x) x(1:end-1),  data.label, 'un', 0));

rpplsWire = cell.empty(size(bndLab,1), 0);
for bund = 1 : size(bndLab,1)
    curBund = cellfun(@(x) regexp(x,bndLab(bund)), {data.label}, 'un', 0);
    curBund = curBund{1};
    curBund = cellfun(@(x) find(x == 1), curBund, 'un', 0);
    curBund = ~cellfun(@isempty, curBund);
    curBund = find(curBund == 1);
    
    
    % rpplWire = zeros(sum(selChan),chanLen*1000);
    % rpplVisu = zeros(sum(selChan),chanLen*1000);
%     rpplsWire{bund,1} = cell.empty(8,0);
    for mw = 1 : size(curBund,1) % LOOP OVER CHANNEL
        
        cfg             = [];
        cfg.channel     = data.label(curBund(mw));         % select MW from current bundle
        curData         = ft_selectdata(cfg, data);
        curData.staging = {ones(1,length(curData.time{1}))};  % add sleep scores
        
        curData.artifacts  = del_sampleinfo(mw)';      % artefacts for the channel of that iteration
        curData.sampleinfo = [1 length(curData.time{1})];
        
        % sets up configuration of ripple detection (saved in variable tfg)
        tfg = rppl_cfg(curData, mw);
        
        % run ripple detection on all wires
        outRipple = Slythm_DetectSpindles_v2(tfg, curData);
        
        rejects = outRipple.falseposRjct.rejects{1};
        numRej  = sum(rejects);
        numRip  = outRipple.evtIndiv.numEvt - numRej;
        
        if numRip >= rppl_exp
            
            rppls = [outRipple.evtIndiv.staTime(~rejects)', outRipple.evtIndiv.endTime(~rejects)'];
            
            rpplsWire{bund,1}{mw,:} = rppls;
                    for rp = 1:size(rppls,1)
%                         rpplWire(mw, rppls(rp,1) : rppls(rp,2)) = 1;
                        rpplVisu(mw, rppls(rp,1) : rppls(rp,1) + 3000) = 1;
                    end
        else
            rpplsWire{bund}{mw,:} = [];
        end
    end
    
    figure('units','normalized','outerposition',[0 0 1 1])
    imagesc(rpplVisu)
    drawnow
    
    % rpplRec = logical(sum(rpplWire,1));
    rpplRec = [];
end 
end % end of function
sub-0002
662 900
S1