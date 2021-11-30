clear
fetch_allSubj;
% 29 + 30
for ia = 1:size(allSubj,1)
    subjID = allSubj{ia};
    bidsID = sub_ID_conversion(subjID, 'yes');
    disp(subjID);
    try
        cd X:/Luca/data
    catch
        cd /media/ldk898/rds-share/Luca/data
    end
    
    mSubject = subjID(1:end-3);
    mSession = subjID(end-1:end);
    
    % if the session name is called 1b then this line prevents an error during cd
    mSubject(regexp(mSubject,'_')) = [];
    if isempty(regexp(mSession,'S', 'ONCE'))
        mSession = ['S', mSession];
    end
    
    cd(mSubject)
    cd(mSession)
    abc = dir; cd(abc(3).name);
    subjPath = cd;
    
    cd([subjPath,'\posDetect'])
    posSpikes = dir('times_*');
    cd('..\negDetect')
    negSpikes = dir('times_*');
    matfiles = [posSpikes; negSpikes];
    allSU = {};
    indexSU = {};
    
    % load electrode location
    cd([subjPath,'\advancedAnalysis\elecLoc']);
    abc = dir('elecLoc_*'); load(abc.name)
    elecLoc = elecLoc(strcmp(elecLoc(:,2),'hipp'),:); % only consider wires inside the hippocampus
    elecLoc(:,1) = cellfun(@(x)regexprep(x,'o_',''),elecLoc(:,1), 'UniformOutput', false);
    
    % load names of cluster that are putative single units (clNames)
    cd([subjPath,'\advancedAnalysis\manualRej']);
    abc = dir('clNames*'); load(abc.name)
    clNames = cellfun(@(x)regexprep(x,'o_',''),clNames, 'UniformOutput', false);

    
    % loop through wires
    for i = 1:size(matfiles,1)
        % extract wire name; e.g.: antHippL1
        wire = matfiles(i).name;
        temp = regexp(wire,'_');
        wire([1:temp(2), end-3:end]) = [];
        
        if ~ismember(wire(1:end-1), elecLoc(:,1)) % is the wire in the hippocampus?
            continue
        end
        
        % positive or negative cluster?
        if ~isempty(regexp(matfiles(i).folder,'pos', 'ONCE'))
            SUsign = 'Pos';
        elseif ~isempty(regexp(matfiles(i).folder,'neg', 'ONCE'))
            SUsign = 'Neg';
        end
        
        load([matfiles(i).folder, filesep, matfiles(i).name], 'spikes', 'cluster_class')
        numCl = max(cluster_class(:,1));
        
        % loop through cluster on that wire
        for ii = 1:numCl
            if ~ismember([wire, SUsign, 'CL',num2str(ii)], clNames) % is the cluster a putative SU and not an artifact?
                continue
            end
            
            spiketimes = cluster_class(cluster_class(:,1)==ii, 2);
            waveshape = spikes(cluster_class(:,1)==ii,:);
            
            % save to global variable
            allSU = [allSU; {wire}, {SUsign}, {spiketimes}, {waveshape}];
        end
    end
    
%     % sanity check
%     if ~isequal(size(allSU,1), size(clNames,1))
%         error('ERROR: Not all good units are included!')
%     else
        save(['X:\Luca\data\allSbj\','allSU_', bidsID,'_', mSession, '.mat'], 'allSU')
%     end
end