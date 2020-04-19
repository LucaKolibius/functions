function rename_clNames(subjID)
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
abc = dir;
cd(abc(3).name)
cd advancedAnalysis
cd allXC
myFile = dir('tableTimestamps*');
load(myFile.name);


% renaming of allXC
for ie = 1: size(tableTimestamps,2)
    extractor = tableTimestamps.Properties.VariableNames{ie};
    if size(extractor,2) ==25
        new_extractor=extractor;
        new_extractor(1:10)='';
        tableTimestamps.Properties.VariableNames{extractor} = new_extractor;
    elseif size(extractor,2) == 21
        new_extractor=extractor;
        new_extractor(1:12)='';
        tableTimestamps.Properties.VariableNames{extractor} = new_extractor;
    elseif size(extractor,2) == 11
        new_extractor=extractor;
        new_extractor(1:2)='';
        tableTimestamps.Properties.VariableNames{extractor} = new_extractor;
    elseif size(extractor,2) == 28
        new_extractor=extractor;
        new_extractor(1:10)='';
        tableTimestamps.Properties.VariableNames{extractor} = new_extractor;
    end
end
save(['tableTimestamps', subjID, '.mat'], 'tableTimestamps');


%% renaming postXCrej
cd ../postXCrej/
myFile = dir ('tableTimestamps_postXC*');
load(myFile.name);
for ie = 1: size(tableTimestamps_postXC,2)
    extractor = tableTimestamps_postXC.Properties.VariableNames{ie};
    if size(extractor,2) ==25
        new_extractor=extractor;
        new_extractor(1:10)='';
        tableTimestamps_postXC.Properties.VariableNames{extractor} = new_extractor;
    elseif size(extractor,2) == 21
        new_extractor=extractor;
        new_extractor(1:12)='';
        tableTimestamps_postXC.Properties.VariableNames{extractor} = new_extractor;
    elseif size(extractor,2) == 11
        new_extractor = extractor;
        new_extractor(1:2) = '';
        tableTimestamps_postXC.Properties.VariableNames{extractor} = new_extractor;
    elseif size(extractor,2) == 28
        new_extractor=extractor;
        new_extractor(1:10)='';
        tableTimestamps_postXC.Properties.VariableNames{extractor} = new_extractor;
    end
end
save(['tableTimestamps_postXC_', subjID, '.mat'],'tableTimestamps_postXC');

%% renaming clNames
try
cd ../manualRej
myFile = dir('clNames*');
load(myFile.name);
if size(clNames,2) < size(clNames,1)
    clNames = clNames';
end

for ie = 1: size(clNames,2)
        extractor = clNames{ie};
        if size(extractor,2) == 25
            new_extractor = extractor;
            new_extractor(1:10) = '';
            clNames(ie) = cellstr(new_extractor);
        elseif size(extractor,2) == 21
            new_extractor = extractor;
            new_extractor(1:12) = '';
            clNames(ie) = cellstr(new_extractor);
        elseif size(extractor,2) == 11
            new_extractor = extractor;
            new_extractor(1:2) = '';
            clNames(ie) = cellstr(new_extractor);
        elseif size(extractor,2) == 28
            new_extractor=extractor;
            new_extractor(1:10)='';
            clNames(ie) = cellstr(new_extractor);
        end
end
save(['clNames',subjID,'.mat'],'clNames');
catch
    disp('No manualRej folder yet');
end
end