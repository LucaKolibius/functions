for ic = 34:size(allSubj,1)
    subjID = allSubj{ic};
    disp(['Looking at Subject: ', subjID])
    if strcmp(subjID(end),'1')
    
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
    cd advancedAnalysis\elecLoc\
    saveCD = cd;

    cd ../..
    abc = dir('CSC_*');
    if isempty(abc)
        abc = dir('Micro_*');
    end
    if isempty(abc)
        warning('raw files not named ''CSC_'' or ''Micro_''!');
    end
    
    % gives me all the different bundle names in that recording
    elecLoc={};
    for ia = 1 : size(abc,1)
        bundleName = abc(ia).name;
        if size(bundleName,2) == 13
            bundleName(1:6) = '';
            bundleName(end-4:end)='';
        else
            bundleName(1:4)='';
            bundleName(end-4:end)='';
        end
        elecLoc{ia,1} = bundleName;
    end
    elecLoc=unique(elecLoc);
    
    % asks for input where each bundle is
    % cortex is the ventral visual stream / other is any other cortical
    % area)
    for ib = 1 : size(elecLoc,1)
        elecLoc{ib,2} = input([ 'Where is the electrode ', elecLoc{ib}, ' (cort / hipp / other)? '], 's');
    end
    
    
    cd(saveCD);
    save(['elecLoc_', subjID(1:end-3), '.mat'], 'elecLoc');
       
    end
end
