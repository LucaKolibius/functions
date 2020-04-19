cd Z:\Luca\data\P08\S1\2017-05-30_17-39-08

par = set_parameters();
par.detection = 'pos';
Get_spikes('all', 'parallel', true, 'par', par);
movefile('*.mat', 'posDetect'); % moves all .mat files that include the just detected spikes into the folder "posDetect"

% neg detection
par.detection = 'neg';
Get_spikes('all', 'parallel', true, 'par', par);
movefile('*.mat', 'negDetect'); % moves all .mat files that include the just detected spikes into the folder "negDetect"

% update pos
cd posDetect/
allDet = dir ('CSC_*_spikes.mat');
for i = 1 : size(allDet, 1)
    load( allDet(i).name, 'index' );
    nameDet = allDet(i).name;

    nameRel = nameDet;
    nameRel(1:4) = []; % delte prefix
    nameRel(end-10:end) = []; % delete suffix
    nameRel = ['times_CSC_', nameRel, '.mat']; % add new prefix & suffix
    load(nameRel);    
    
    cluster_class(:,2) = index; % update with correct timestamps
    save(nameRel , 'spikes', 'cluster_class', 'par', 'gui_status', 'Temp', 'forced', 'inspk', 'ipermut');
    
    clearvars -except allDet
end

% update neg
cd ../negDetect/
allDet = dir ('CSC_*_spikes.mat');
for i = 1 : size(allDet, 1)
    load( allDet(i).name, 'index' );
    nameDet = allDet(i).name;

    nameRel = nameDet;
    nameRel(1:4) = []; % delte prefix
    nameRel(end-10:end) = []; % delete suffix
    nameRel = ['times_CSC_', nameRel, '.mat']; % add new prefix & suffix
    load(nameRel);    
    
    cluster_class(:,2) = index; % update with correct timestamps
    save(nameRel , 'spikes', 'cluster_class', 'par', 'gui_status', 'Temp', 'forced', 'inspk', 'ipermut');
    
    clearvars -except allDet
end
