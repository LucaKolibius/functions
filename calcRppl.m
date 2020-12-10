function [rpplsWire, rpplLen, bndLab] = calcRppl (data)

chanLen   = cell2mat(cellfun(@length, [data.time], 'un', 0)) / 1000; % timepoint length divided by sampling freq gives channellength in seconds
% rppl_exp  =  chanLen * (0.023 + 0.00545);           % number of ripples in white noise of same length
rppl_exp = 0;

% NAMES OF THE BUNDLES
bndLab = unique(cellfun(@(x) x(1:end-1),  data.label, 'un', 0));



% rpplsWire = cell.empty(size(data.trial,2), 0);
rpplLen   = cell(1,size(data.trial,2));
rpplsWire = cell(1,size(data.trial,2));
for trl = 1 : length(rpplsWire)
    rpplsWire{trl} = zeros(1,size(data.trial{trl},2)); % gets overwritten later for no ripples
end

for bund = 1 : size(bndLab,1) % currently 1!
    curBund = cellfun(@(x) regexp(x,bndLab(bund)), {data.label}, 'un', 0);
    curBund = curBund{1};
    curBund = cellfun(@(x) find(x == 1), curBund, 'un', 0);
    curBund = ~cellfun(@isempty, curBund);
    curBund = find(curBund == 1);
    
    for mw = 1 : size(curBund,1) % LOOP OVER CHANNEL | currently 1!
        
        cfg             = [];
        cfg.channel     = data.label(curBund(mw));         % select MW from current bundle
        curData         = ft_selectdata(cfg, data);
        
        %% DEFINE ARTEFACTS
        for trl = 1:size(curData.trial,2)
            curData.staging{1,trl}           = ones(1, round(chanLen(trl)*1000));  % add sleep scores
            [~, curData.artifacts{1,trl}]    = iqrAR(curData.trial{trl},0);
            %                 curData.sampleinfo(trl,:)        = [1 chanLen(trl)*1000];
        end
        
        % sets up configuration of ripple detection (saved in variable tfg)
        tfg = rppl_cfg(curData, mw);
        
        % run ripple detection on all wires
        outRipple = Slythm_DetectSpindles_v2(tfg, curData);
        
        for trl = 1:size(curData.trial,2)
            
            if outRipple.evtSummary(trl).numEvt == 0 % if there are no ripples, rejection was skipped in the ripple algorithm
                numRip(trl) = 0;
                
            else % most cases
                rejects = outRipple.falseposRjct.rejects{trl};
                numRej  = sum(rejects);
                numRip(trl)  = outRipple.evtIndiv(trl).numEvt - numRej;
            end
            
            if numRip(trl) > 0
                
                switch isempty(rejects)
                    case 0
                        rppls = [outRipple.evtIndiv(trl).staTime(~rejects)', outRipple.evtIndiv(trl).endTime(~rejects)'];
                    case 1 % no rejects
                        rppls = [outRipple.evtIndiv(trl).staTime', outRipple.evtIndiv(trl).endTime'];
                end
                
                rpplLen{trl} = ( rppls(:,2) - rppls(:,1) + ones(size(rppls,1),1) )';
                
                %                 rpplsWire{bund,1}{mw,:} = rppls;
                for rp = 1:size(rppls,1)
                    rpplsWire{trl}(rppls(rp,1) : rppls(rp,2)) = 1;
                end
            else
                rpplsWire{trl} = []; % no ripples in this trial
            end
            
            
        end % END OF TRIAL
        
        
    end % END OF MW
end % END OF BUNDLE
end % END OF FUNCTION