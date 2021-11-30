clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\ESN_code_data\code_ESN\data\allSpks.mat', 'allSpks')

for spk = 1 : length(allSpks)
    
    subj = sub_ID_conversion(allSpks(spk).bidsID, 'yes');
    sesh = allSpks(spk).sesh;
    
    if strcmp(subj, 'P7_ERL')
        subj = 'P07ERL';
    end
    
    try
        cd(['\\analyse4.psy.gla.ac.uk\project0309\Luca\data', filesep, subj, filesep, sesh])
    catch
        cd(['/analyse/Project0309/Luca/data', filesep, subj, filesep, sesh])
    end
    
    abc = dir; cd(abc(3).name);
    p2d = [cd, filesep];
    
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, himiDx] = loadLogs(p2d, 1);
   
    allSpks(spk).himiDx = himiDx;
end
save('\\analyse4.psy.gla.ac.uk\project0309\Luca\ESN_code_data\code_ESN\data\allSpks.mat', 'allSpks')

