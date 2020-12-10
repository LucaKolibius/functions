function [lib] = checkServerData_sub(lib, fieldnam)
disp(cd);
allCont = dir(cd);
allDat = allCont(~[allCont.isdir]);
allFold = allCont([allCont.isdir]);

allFold(strcmp({allFold.name}, 'code4Mem')) = []; % skip code4Mem folder (only experimental code in here)

%% ADD DATA TO LIBRARY
for dat = 1 : size(allDat,1)
    posi = size(lib,2)+1;
    for f = 1:6
        lib(posi).(fieldnam{f}) = allDat(dat).(fieldnam{f});
    end
end

for folder = 1 : size(allFold,1)
    if strcmp(allFold(folder).name, '.') || strcmp(allFold(folder).name, '..')
        continue
    end
    
%     changedDir = 0;
%     while changedDir == 0
%     try
        cd([allFold(folder).folder, filesep, allFold(folder).name])
%         changedDir = 1;
%     catch
%     end
%     end
    
    [lib] = checkServerData_sub(lib, fieldnam);
end

end