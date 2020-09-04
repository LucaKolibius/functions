function [resLen, resNum] = spkRppl_anal_sub1(allSU, trigALL, trigIU, tw)

for suIt = 1 : size(allsU,1)
    disp(suIt)
    
    % GET CURRENT bidsID AND sesh
    bidsID = allSU(suIt).bidsID;
    sesh   = allSU(suIt).sesh;
    su = allSU(suIt).su;
    wirename = trigALL(suIt).wirename;
    bundlename = wirename(1:end-1);
    
    % REPEATING WIRES
    if su > 1
        if      and(strcmp( trigALL(su-1).bidsID,            trigALL(su).bidsID            ), ...
                and(strcmp( trigALL(su-1).sesh,              trigALL(su).sesh              ), ...
                    strcmp( trigALL(su-1).wirename(1:end-1), trigALL(su).wirename(1:end-1) )))
            
            continue
            
        end
    end
   
    % FIND OUT WHICH TRIALS ARE INDEXED IN THAT BUNDLE
    allBund = {trigALL.wirename};
    allBund = cellfun(@(x) x(1:end-1), allBund, 'un', 0)';
    
   labs   = trigALL.wirename;

    
    lfpNam = lfpNam - noTrig - noHits;                          % index all the names of LFPs
    LFPlab = labs(logical(lfpNam));                             % names of all the LFPs (MB1, MB2, MB3...)
    bndLab = unique(cellfun(@(x) x(1:end-1), LFPlab, 'un', 0)); % only the unique bundles
    
    for bnd = 1 : size(bndLab,1)
        
        % which rows of lfp data do we need from data_micro?
        microIdx = cellfun(@(x) regexp(x,bndLab(bnd)), {labs}, 'un', 0);
        microIdx = microIdx{1};
        microIdx = cellfun(@(x) find(x == 1), microIdx, 'un', 0);
        microIdx = ~cellfun(@isempty, microIdx);
        
        % get the LFP from data_micro
        dat = data_micro.trial{1}(microIdx,:);
        
        % FIND THE TRIALS OF THAT BUNDLE THAT ARE INDEXING (IF ANY)
        % extract all SU that fire during that session + bundle
        allBund  = {trigALL.wirename};
        allBund  = cellfun(@(x) x(1:end-1), allBund, 'un', 0);
        sameBund = and(and(contains({trigALL.bidsID}, bidsID), strcmp({trigALL.sesh}, sesh) ), contains(allBund, bndLab(bnd))); % same subject + session + bundle
        
        idxTrl = any(vertcat(trigALL(sameBund).idxTrl),1); % these are the trials that are indexed in that wire
        




end