baseDir = '\\analyse4.psy.gla.ac.uk\project0309\Luca\data\';
for su = 1:length(allSpks)
    bidsID = allSpks(su).bidsID;
    sesh   = allSpks(su).sesh;
    subjID = sub_ID_conversion(bidsID, 1);
    
    cd([baseDir,subjID, filesep, sesh])
    abc = dir('20*'); cd(abc.name);
    
    [~, ~, ~, ~, ~, retTrigger, ~, retRT, ~, ~] = loadLogs([cd, filesep], 1);
    
    allSpks(su).retRT = retRT;
end
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ.mat', 'allSpks');