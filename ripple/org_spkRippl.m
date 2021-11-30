% go through tfg settings with a fine comb

% artDat = artDat( cellfun(@(x) isempty(x), regexp({artDat.name}, '_th-8')) ); % without files that have the '_th-8' extension 

%% sub 1003_S1 makes problems with indexing
clear; clc;
IU = struct; numIU = 0;
SU = struct; numSU = 0;
sanCheck      = 0; % checks how many units are skipped. has to be equal to 13 / the number of IU
suNum         = 0;
spksRppl_norm = [];
server        = 0; % 0 for local, 1 for server

switch server
    case 0
        %         lfpFold = 'X:\George\Analysis\Data Continuous';                                                    % LFP
        lfpFold = 'X:\Luca\data\microLFP';
        artFold = 'Z:\hanslmas-ieeg-compute\George\Analysis\Artefact Rejection\Data Continuous 1000Hz';    % ARTEFACTS
        suDir   = dir('X:\Luca\data\allSbj\allSU_sub-*.mat');                                              % ALL SU
        load('X:\Luca\indexSUvisu\varWindow_ENC_-1-113_RET_-1113_DP_p99_th1\allIU.mat', 'allIU');          % WHICH ONES ARE IU
    case 1
        lfpFold = '/media/ldk898/rds-share/Luca/server/George';                                            % LFP
        artFold = '/media/ldk898/rds-share/Luca/server_ripples/Artefact Rejection/Data Continuous 1000Hz'; % ARTEFACTS
        suDir   = dir('/media/ldk898/rds-share/Luca/server/data/allSbj');                                  % ALL SU
        load('/media/ldk898/rds-share/Luca/server_ripples/allIU.mat', 'allIU');                            % WHICH ONES ARE IU
end

for sbj = 1:size(suDir,1)      % LOOP OVER SESSIONS
    
    switch server
        case 0
            load([suDir(1).folder, filesep, suDir(sbj).name], 'allSU'); % LOAD ALL SPIKES OF THAT SESSION
        case 1
            
    end
    
    bidsID = suDir(sbj).name; bidsID = bidsID(7:14);
    sesh   = suDir(sbj).name; sesh   = sesh(16:18) ; sesh = regexprep(sesh,'\.','');

    
    for su = 1 : size(allSU,1) % LOOP OVER SU
        wireName = allSU{su,1};
        su_spks  = round(allSU{su, 3}/32); % spiketimes of that SU

        %% CHECK IF SU IS A IU
        if any( and( and( contains(allIU(:,1), bidsID), contains(allIU(:,2), sesh) ), ([allIU{:,3}] == su)') )
            sanCheck = sanCheck + 1;
            isIU     = 1;
            numIU    = numIU + 1;
        else
            isIU     = 0;
            numSU    = numSU + 1;
        end

        %% get the lfp data
        try
%             load([lfpFold, filesep, bidsID, '_', sesh, '_micro_RAW_DS-1000_SPK-INT.mat'], 'data_micro')
            load([lfpFold, filesep, bidsID, '_', sesh, '_onlyMicroLFP_RAW_1000DS_SPKINT.mat'], 'data')
        catch
            skipped(sbj) = 1;
            continue
        end

        suNum = suNum +1;


        %% run calculations and output rpplRec
        if su == 1
            [rpplRec, rpplWire] = calcRppl (data, bidsID, sesh, artFold);
        elseif su > 1
            if ~strcmp(wireName(1:end-1), allSU{su-1,1}(1:end-1))
                [rpplRec, rpplWire] = calcRppl (data_micro, bidsID, sesh, artFold);
            end
        end
        
        %% save calculations in variable
        switch isIU
            case 1
                IU(numIU).rpplRec  = rpplRec;
                IU(numIU).spks     = su_spks;
                IU(numIU).bidsID   = bidsID;
                IU(numIU).sesh     = sesh;
                IU(numIU).su       = su;
                IU(numIU).rpplWire = rpplWire;
            case 0
                SU(numSU).rpplRec  = rpplRec;
                SU(numSU).spks     = su_spks;
                SU(numSU).bidsID   = bidsID;
                SU(numSU).sesh     = sesh;
                SU(numSU).su       = su;
                SU(numSU).rpplWire = rpplWire;
        end
    end
    
    % save variable
%     switch server
%         case 0
%             save('X:\Luca\ripple and IU\rpplRec', 'IU', 'SU', '-v7.3');
%         case 1
%             save('/media/ldk898/rds-share/Luca/server_ripples/rpplRec', 'IU', 'SU');
%     end
end
