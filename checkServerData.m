clear
castlesAr = 'X:\Archive';
rdsAr     = 'Z:\Archive';
trueAr    = 'Y:\Archive';

cd \\analyse4.psy.gla.ac.uk\project0309\Luca
%% GET A LIBRARY OF ARCHIVE CONTENTS
for it = 1 %: 3
    
%     switch it
%         case 1
%             cd(castlesAr)
%         case 2
%             cd(rdsAr)
%         case 3
%             cd(trueAr)
%     end
    
    lib    = struct('name', [], 'folder', [], 'date', [], 'bytes', [], 'isdir', [], 'datenum', []); % preallocate library
    lib(1) = [];
    fieldnam   = {'name'; 'folder'; 'date'; 'bytes'; 'isdir'; 'datenum'};
    
    [lib] = checkServerData_sub(lib, fieldnam);
    
    switch it
        case 1
            castlesLib = lib;
        case 2
            rdsLib = lib;
        case 3
            archLib = lib;
    end
end

%% compare castles and archive
allArchName  = {archLib.name};
allArchBytes = [archLib.bytes];
isOnArchive  = zeros(1, size(castlesLib,2));
for dat = 1 : size(castlesLib,2)
    disp(dat)
    name = castlesLib(dat).name;
    bytes = castlesLib(dat).bytes;
    
    isOnArchive(dat) = any(strcmp(name, allArchName) & bytes == allArchBytes);
    
end

sum(~isOnArchive)
missing = find(isOnArchive == 0);
% missing(4)

missingDat = castlesLib(missing);

checkNam = cellfun(@(x) regexp(x, 'P02bak'), {missingDat.folder}, 'un', 0);
checkNam = cellfun(@isempty, checkNam, 'un', 0); checkNam = cell2mat(checkNam);
sum(checkNam)

checkNam = cellfun(@(x) regexp(x, 'P02'), {archLib.folder}, 'un', 0);
checkNam = cellfun(@isempty, checkNam, 'un', 0); checkNam = cell2mat(checkNam); checkNam = logical(~checkNam);
findinArch = archLib(checkNam);

    
