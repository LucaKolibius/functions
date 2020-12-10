clear
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_fixed2.mat')
addpath(genpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\functions'));
addpath(genpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\toolboxes\wave_clus-master'));
addpath(genpath('\\analyse4.psy.gla.ac.uk\project0309\Luca\toolboxes\TREBER'))
% checkchannelLDK(1)
prevWire = [];
for su = [583, 586, 590, 595, 619]%1:length(allSpks)
    disp(su)
    newWire = [allSpks(su).bidsID, allSpks(su).sesh, allSpks(su).wirename, num2str(allSpks(su).isPos)];
    
     if strcmp(prevWire, newWire)
        continue
     else
         prevWire = newWire;
     end
   
      
%     disp(newWire);
    cd('\\analyse4.psy.gla.ac.uk\project0309\Luca\data')
    bidsID = allSpks(su).bidsID;
    sesh   = allSpks(su).sesh;
    sr     = allSpks(su).initSR;
    wirename = allSpks(su).wirename;
    bundlename = allSpks(su).bundlename;
    isPos = allSpks(su).isPos;
    
    switch isPos
        case 0
            posNeg = 'neg';
        case 1
            posNeg = 'pos';
    end
    
    sbj = sub_ID_conversion(bidsID,1);
    cd(sbj); cd(sesh); abc = dir('2*'); cd(abc.name);
    
    if exist('detect_newSR', 'file') == 7
        cd('detect_newSR');
    else
        switch allSpks(su).isPos
            case 0
                cd('negDetect');
            case 1
                cd('posDetect');
        end
    end
    
    det = dir(['*', wirename, '*_spikes.mat']); load(det.name, 'par'); parDet = par.sr;
    clus =  dir(['times_*', wirename, '*.mat']); load(clus.name, 'par', 'cluster_class'); parClus = par.sr;
    clear par;
    
        loadNumSpks = sum(cluster_class(:,1)==allSpks(su).su);

        %% CHECKING FOR CORRECT NUMBER OF SPIKES (THIS WAS FAULTY EARLIER BECAUSE I DID NOT SEPARATE NEG AND POS IN detect_newSR
        if loadNumSpks~=length(allSpks(su).spks)
            cd('..\detect_newSR2\')
            clus =  dir(['times_*', wirename, '_', posNeg,'*.mat']); load(clus.name, 'par', 'cluster_class');
            loadNumSpks = sum(cluster_class(:,1)==allSpks(su).su);
        end
        
        if loadNumSpks~=length(allSpks(su).spks)
            disp(loadNumSpks);
            disp(length(allSpks(su).spks));
        end
        
    %% CHECKING FOR SAMPLING RATE
%     if ~isequal(parDet, parClus, sr) % is there a difference?
%         cd('..') % go to raw data
%         switch isequal(parDet, sr) % does the detection-sr match the real-sr?
%             case 0 % redo detection + clustering
%                 disp('REDOING DETECTION AND CLUSTERING!')
%                 
%                 %% REDO DETECTION
%                 par = set_parameters();
%                 par.detection = num2str(allSpks(su).isPos);
%                 par.detection = regexprep(par.detection, '0', 'neg');
%                 par.detection = regexprep(par.detection, '1', 'pos');
%                 par.sr = sr;
%                 thisWire = dir(['*',wirename, '*.ncs']);  
%                 if size(thisWire,1) > 1; error('too many files'); end
%                 
%                 Get_spikes(thisWire.name, 'parallel', false, 'par', par);
%                 
%                 %% RENAME FILE & MOVE TO detect_newSR2
%                 thisWire = dir(['*',wirename, '*_spikes.mat']);
%                 thisWire = thisWire.name;
%                 
%                 switch allSpks(su).isPos
%                     case 0
%                         newname = regexprep(thisWire, '_spikes.mat', '_neg_spikes.mat');
%                     case 1
%                         newname = regexprep(thisWire, '_spikes.mat', '_pos_spikes.mat');
%                 end
%                 movefile(thisWire, newname)
%                 mkdir('detect_newSR2')
%                 movefile(newname, 'detect_newSR2'); % moves all .mat files tha
%                 
%                 %% REDO CLUSTERING
%                 cd('detect_newSR');
%                 thisWire = dir(['*',wirename, '_',posNeg, '_spikes.mat']);
%                 if size(thisWire,1) > 1; error('too many files'); end
%                 
%                 Do_clustering(thisWire.name, 'parallel', false, 'make_plots', false);
%                 
%                 moveTo = cd;
%                 
%             case 1 % used the correct SR during detection, only redoing clustering                        
%                          %% GO TO THE RIGHT POLARITY FOLDER
%                 switch allSpks(su).isPos
%                     case 0
%                         cd('negDetect');
%                     case 1
%                         cd('posDetect');
%                 end
%     
%                 thisWire = dir(['*',wirename, '_',posNeg, '_spikes.mat']);
%                 if size(thisWire,1) > 1; error('too many files'); end
%            
%                 Do_clustering([thisWire.folder, filesep, thisWire.name], 'parallel', false, 'make_plots', false); % already inherits correct par variable
%                 moveTo = cd;
%                 moveTo(end-8:end) = [];
%                 moveTo = [moveTo, 'detect_newSR'];
%                 movefile(['times_*', wirename, '.mat'], moveTo); % moves all .mat files tha
%         end
%                 % recluster
%                 sprintf('CL1 #%.f || CL2 #%.f || CL3 #%.f', sum(cluster_class(:,1)==1), sum(cluster_class(:,1)==2), sum(cluster_class(:,1)==3))
%                 cd(moveTo);
%                 abc = dir(['times_*', wirename, '.mat']); !! posNeg
%                 if size(abc,1) > 1; error('too many files'); end
%                 
%                 figure(1)
%                 mhandle = wave_clus(abc.name); % loads results from automatic clustering
%                 set(mhandle,'WindowStyle','normal'); % undock
%                 %     set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
%                 set(mhandle,'units','normalized','OuterPosition',[0.9958   -0.0815    1.0083    1.0889]) % fullscreen
%                 close(1);
%                 
%                 nxt = input('Load next SU (y)? ', 's');
%                 while nxt~='y'
%                     nxt = input('How about now: Load next SU (y)? ', 's');
%                 end
%                 allSpks(su).fixed = 1;
%                 
%                 load(abc.name, 'cluster_class')


% %%
%                 newCl = max(cluster_class(:,1));
%                 for cl = 1:newCl
%                     allSpks(size(allSpks,2)+1).bidsID     = bidsID;
%                     allSpks(size(allSpks,2)).sesh       = sesh;
%                     allSpks(size(allSpks,2)).wirename   = wirename;
%                     allSpks(size(allSpks,2)).bundlename = bundlename;
%                     allSpks(size(allSpks,2)).su         = cl;
%                     allSpks(size(allSpks,2)).isPos      = isPos;
%                     allSpks(size(allSpks,2)).encTrigger = allSpks(su).encTrigger;
%                     allSpks(size(allSpks,2)).retTrigger = allSpks(su).retTrigger;
%                     allSpks(size(allSpks,2)).hitsIdx    = allSpks(su).hitsIdx;
%                     
%                     spks = cluster_class(cluster_class(:,1)==cl, 2);
%                     spks = spks/sr*1000;
%                     allSpks(size(allSpks,2)).spks       = spks;
%                     allSpks(size(allSpks,2)).iu         = [];
%                     allSpks(size(allSpks,2)).idxTrl     = [];
%                     allSpks(size(allSpks,2)).initSR     = sr;
%                     allSpks(size(allSpks,2)).fixed      = 2;
%                 end
% 
% 
% 
%                 save('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allSpksHZ_fixed', 'allSpks');
%     end
end

 