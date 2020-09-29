function [rpplRec, rpplsWire, bndLab] = calcRppl (data, bidsID, sesh, artFold, ortho)

chanLen   =  length( data.time{1} ) ./ 1000;  % timepoint length divided by sampling freq 1024 Hz gives channellength in seconds
rppl_exp  =  chanLen * (0.023 + 0.00545);           % number of ripples in white noise of same length

% reject artefact wires
load([artFold, filesep, bidsID, '_', regexprep(sesh, 'S1b', 'S1'), '_micro_th-8.mat'], 'del_sampleinfo', 'delIx')

% NAMES OF THE BUNDLES
bndLab = unique(cellfun(@(x) x(1:end-1),  data.label, 'un', 0));

% ORTHOGONALIZE 
if ortho == 1
    [data, success] = orthogonVec(data);
end

rpplsWire = cell.empty(size(bndLab,1), 0);
for bund = 1 : size(bndLab,1)
    curBund = cellfun(@(x) regexp(x,bndLab(bund)), {data.label}, 'un', 0);
    curBund = curBund{1};
    curBund = cellfun(@(x) find(x == 1), curBund, 'un', 0);
    curBund = ~cellfun(@isempty, curBund);
    curBund = find(curBund == 1);
        
    
    % rpplWire = zeros(sum(selChan),chanLen*1000);
    % rpplVisu = zeros(sum(selChan),chanLen*1000);
    % rpplsWire{bund,1} = cell.empty(8,0);
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
        
        if outRipple.evtSummary.numEvt == 0 % if there are no ripples, rejection was skipped in the ripple algorithm
            numRip = 0;
            
        else % most cases
            rejects = outRipple.falseposRjct.rejects{1};
            numRej  = sum(rejects);
            numRip  = outRipple.evtIndiv.numEvt - numRej;
        end
        
        if numRip >= rppl_exp
            
            rppls = [outRipple.evtIndiv.staTime(~rejects)', outRipple.evtIndiv.endTime(~rejects)'];
            
            rpplsWire{bund,1}{mw,:} = rppls;
            for rp = 1:size(rppls,1)
                rpplWire(mw, rppls(rp,1) : rppls(rp,2)) = 1;
%                 rpplVisu(mw, rppls(rp,1) : rppls(rp,1) + 3000) = 1;
            end
        else
            rpplsWire{bund}{mw,:} = [];
        end
    end
    
%     figure('units','normalized','outerposition',[0 0 1 1])
%     imagesc(rpplVisu)
%     drawnow
    
    % rpplRec = logical(sum(rpplWire,1));
    rpplRec = []; % no longer in use (this was taking all the ripples on the bundle into account; independent on mw)
end 
end % end of function
% sub-0002
% 662 900
% S1