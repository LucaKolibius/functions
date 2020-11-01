function xc_postClust(subjID)
mSubject = subjID(1:end-3);
mSession = subjID(end-1:end);

try
    cd X:/Luca/data
catch
    %     cd /media/ldk898/rds-share/Luca/data
    cd /castles/nr/projects/h/hanslmas-ieeg-compute\Luca\data
    
end

cd(mSubject)
cd(mSession)
abc = dir;
cd(abc(3).name)
mkdir('advancedAnalysis');

% pos spikes
% create a table with the spike timestamps for all wires, subdivided into
% cluster (includes 0 cluster with artifacts)
cd posDetect
allSpks = dir('times_*'); % load all wires
tableTimestamps=table.empty;
for ia = 1:size(allSpks,1)
    load(allSpks(ia).name, 'cluster_class'); % load one wire
    
    wirename = allSpks(ia).name;
    wirename(end-3:end)='';
    if size(wirename,2)==8 %% P01ERL has a different raw datafile name
        wirename(1:6)='';
    elseif size(wirename,2)==19 || size(wirename,2)==20
        wirename(1:10)='';
    elseif size(wirename,2)==15 % for P07ERL (wirename example: 'times_Micro_MA1')
        wirename(1:12)=[];
    end
    wirename(end+1:end+3)='Pos';
    
    allClust = max(cluster_class(:,1)); % how many cluster are in one wire
    
    % create table variable-names
    varname = cell.empty;
    for ib = 0:allClust
        varname(1,size(varname,2)+1) = cellstr([wirename, sprintf('CL%d', ib)]);
    end
    
    % create empty table with variable-names
    tempTable = cell2table(cell(1,size(varname,2)), 'VariableNames',varname);
    
    % insert spiketimes into table
    for ic = 0:allClust
        tempTable{1 , ic+1} = {cluster_class(cluster_class(:,1)==ic, 2)}; % ic+1 because it starts at "0"
    end
    
    tableTimestamps = [tableTimestamps tempTable];
    clear tempTable
end

%% negative spikes
cd ../negDetect
allSpks = dir('times_*'); % load all wires

% create a table with the spike timestamps for all wires, subdivided into
% cluster (includes 0 cluster with artifacts)
for ia = 1:size(allSpks,1)
    load(allSpks(ia).name, 'cluster_class'); % load one wire
    
    wirename = allSpks(ia).name;
    wirename(end-3:end)='';
    if size(wirename,2)==8 %% P01ERL has a different raw datafile name
        wirename(1:6)='';
    elseif size(wirename,2)==19 || size(wirename,2)==20
        wirename(1:10)='';
    end
    wirename(end+1:end+3)='Neg';
    
    allClust = max(cluster_class(:,1)); % how many cluster are in one wire
    
    % create table variable-names
    varname = cell.empty;
    for ib = 0:allClust
        varname(1,size(varname,2)+1) = cellstr([wirename, sprintf('CL%d', ib)]);
    end
    
    % create empty table with variable-names
    tempTable = cell2table(cell(1,size(varname,2)), 'VariableNames',varname);
    
    % insert spiketimes into table
    for ic = 0:allClust
        tempTable{1 , ic+1} = {cluster_class(cluster_class(:,1)==ic, 2)}; % ic+1 because it starts at "0"
    end
    
    tableTimestamps = [tableTimestamps tempTable];
    clear tempTable
end

%% find xc between all clusters
% make histogram so the delay between two cells is fixed
BW = 3.2*7; % 0.7ms bin width
startHist = 0; % start histogram
tempMax = [];
for ia = 1:size(tableTimestamps,2)
    if isempty(tableTimestamps{1,ia}{1})
        continue
    end
    tempMax(1,ia)=max(tableTimestamps{1,ia}{1});
end
endHist = []; % end histogram
endHist = max(tempMax); % find last recorded spike

dt = [startHist: BW :endHist]; % 0.1ms bins / linear space for histogram
maxlag = 3; % 2.1ms

% preallocate deleteClust
deleteClust = zeros(7, size(tableTimestamps,2) * size(tableTimestamps,2));
indexCount = 0;

for ia=1:size(tableTimestamps,2)
    a1=hist(tableTimestamps{1,ia}{1},dt);
    for ib=1:size(tableTimestamps,2)
        if ia==ib || ia>ib % im not interested in the autocorrelation and the matrix is symetrical
            continue
        end
        
        indexCount = indexCount+1;
        
        % cross-correlation between two noise cluster are not interesting
        varname1 = tableTimestamps.Properties.VariableNames{ia};
        varname2 = tableTimestamps.Properties.VariableNames{ib};
        if strcmp(varname1(end-2:end),'CL0') && strcmp(varname2(end-2:end),'CL0')
            continue
        end
        
        disp([num2str(ia), ' - ', num2str(ib)]);
        a2=hist(tableTimestamps{1,ib}{1},dt);
        
        if sum(a1)==0 || sum(a2)==0 % one wire had no CL0 entries. the xc produced an error
            continue
        end
        
        [xc,indx] = xcorr(a1, a2, maxlag, 'coeff'); % cross correlation
        
        % indexing clusters that exceed a specific correlation
        deleteClust(1,indexCount) = ia; % first cluster
        deleteClust(2,indexCount)   = ib; % second cluster
        deleteClust(3,indexCount)   = max(xc); % record highest cross-correlation
        temp = indx(xc==max(xc)); % record index where cross-correlation is maximal
        temp2 = temp(randi(size(temp,2))); % if I have two values for the peak position, choose one at random (will likely only happen for very small XC)
        deleteClust(4,indexCount)   = temp2;
    end
end

% trim columns in which the first row is a 0 (comes from preallocation in L103)
deleteClust(:,deleteClust(1,:)==0)=[];

% save variables
cd ../advancedAnalysis/
mkdir('allXC');
cd allXC
save(['deleteClust',subjID ,'.mat'], 'deleteClust'); % has all the xc between all clusters
save(['tableTimestamps', subjID, '.mat'], 'tableTimestamps'); % 1xn with n being the number of cluster. each cell has all spiketimes in samples

% only work with the xc above th
th = 0.20;
mindx = deleteClust(3,:)>=th;
deleteClust_th = deleteClust(:,mindx);

% save new deleteClust under a new variable name
save(['deleteClust_th',subjID,'.mat'],'deleteClust_th'); % thresholded (only crosscorrelations abova th
mCD = cd; % the allXC folder
end
