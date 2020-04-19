% outdated
function animalCueSimil
clear
load('X:\Luca\anCueNames.mat', 'anCueNames');
nperm = 100000;


simil = zeros(1,nperm);
for ia=1:nperm
    idx = randperm(size(anCueNames,1),2);
    [vec1, vec2] = anCueNames{idx,3};
    simil(ia) = cosSimil(vec1,vec2);
end
% hist(simil)
prctile(simil,90)

%%
datapath = 'X:\Luca\indexSUvisu\DP_p95_th2';
cd(datapath)
allSU = dir('SU*');
allEmpSimi = [];
for indexNeuron = 1:size(allSU,1)
    cd(datapath)
    load(allSU(indexNeuron).name)
    
    % cd
    try
        cd X:/Luca_old/data
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
    
    % logfile
    p2d = cd;
    p2d(end+1)='\';
    
    ffTrls = find(ffTrls==1);
    [~, ~, ~, ~, ~, ~, ~, ~, animalCuesFF] = loadLogs(p2d, ffTrls);
    ffAnimalVec = anCueNames(ismember(anCueNames(:,1), animalCuesFF),3);
    
    ppTrls = find(ppTrls==1);
    [~, ~, ~, ~, ~, ~, ~, ~, animalCuesPP] = loadLogs(p2d, ppTrls);
    ppAnimalVec = anCueNames(ismember(anCueNames(:,1), animalCuesPP),3);
    
    %% this part extracts the minimal similarity between two animal cues in case the SU is sensitive towards multiple 
    empSimi = [];
    for ff = 1:length(ffTrls)
        for pp = 1:length(ppTrls)
            empSimi = [empSimi, cosSimil(ffAnimalVec{ff}, ppAnimalVec{pp})];
        end
    end
    empSimi = min(empSimi);
    allEmpSimi = [allEmpSimi empSimi];
    
    %% now extract the trials during which an animal with closer similarity was shown during encoding and retrieval
    % load all shown animal cues from logfile (create new function)
    
    % find animal cues that are closer to the index animal cues (which ones to choose if there are multiple?)
end
end