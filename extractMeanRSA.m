RSAmean_h_all = [];
RSAmean_c_all = [];
RSAmean_o_all = [];
fetch_allSubj
for i = 1 : size(allSubj,1)
    subjID = allSubj{i};
    clear RSAmean_h RSAmean_c RSAmean_o
    
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
    cd advancedAnalysis\
    cd RSA
    abc = dir('RSAmean_*');
    load(abc.name)
    
    %%
    if ~isempty(RSAmean_h) % in case there are not enough hippocampal single units
        if size(RSAmean_h,1) == 3 % in case there are not enough misses
            RSAmean_h(4:6, 1:8) = nan;
            RSAmean_h(4:6, 9  ) = [1;2;3];
            RSAmean_h(4:6, 10 ) = [0;0;0]; % redundant
        end
        RSAmean_h(1:size(RSAmean_h,1),end+1) = i;
        RSAmean_h_all = [RSAmean_h_all; RSAmean_h];
    end
    
    if ~isempty(RSAmean_c)
        if size(RSAmean_c,2) == 3
            RSAmean_c(4:6, 1:8) = nan;
            RSAmean_c(4:6, 9  ) = [1;2;3];
            RSAmean_c(4:6, 10 ) = [0;0;0]; % redundant
        end
        RSAmean_c(1:size(RSAmean_c,1),end+1) = i;
        RSAmean_c_all = [RSAmean_c_all; RSAmean_c];
    end
    
    if ~isempty(RSAmean_o)
        if size(RSAmean_o,2) == 3
            RSAmean_o(4:6, 1:8) = nan;
            RSAmean_o(4:6, 9  ) = [1;2;3];
            RSAmean_o(4:6, 10 ) = [0;0;0]; % redundant
        end
        RSAmean_o(1:size(RSAmean_o,1),end+1) = i;
        RSAmean_o_all = [RSAmean_o_all; RSAmean_o];
    end
end

% clearvars -except RSAmean_h_all RSAmean_c_all RSAmean_o_all


% RSAmean_h_1Hi = RSAmean_h_all(1:6:end,:);
% RSAmean_h_2Hi = RSAmean_h_all(2:6:end,:);
% RSAmean_h_3Hi = RSAmean_h_all(3:6:end,:);
% RSAmean_h_1Mi = RSAmean_h_all(4:6:end,:);
% RSAmean_h_2Mi = RSAmean_h_all(5:6:end,:);
% RSAmean_h_3Mi = RSAmean_h_all(6:6:end,:);
% idx(:,1)=~isnan(RSAmean_h_1Mi(:,1));
% idx(:,2)=~isnan(RSAmean_h_2Mi(:,1));
% idx(:,3)=~isnan(RSAmean_h_3Mi(:,1));
% idx=find(sum(idx,2)==3);
% 
% newRSAmeanH = [RSAmean_h_1Hi(:,1:8) RSAmean_h_1Mi(:,1:8) RSAmean_h_2Hi(:,1:8) RSAmean_h_2Mi(:,1:8) RSAmean_h_3Hi(:,1:8) RSAmean_h_3Mi(:,1:8)];

% reformat hippocampus
x(:,1) = (1:6:246)'; % WT hit
x(:,2) = (4:6:246)'; % WT miss
x(:,3) = (2:6:246)'; % WC hit
x(:,4) = (5:6:246)'; % WC miss
x(:,5) = (3:6:246)'; % BC hit
x(:,6) = (6:6:246)'; % BC miss

newRSAmeanH = [];
for ix = 1:8
    for ib = 1:6
    temp1 = RSAmean_h_all(:,ix);
    newRSAmeanH(:,size(newRSAmeanH,2)+1) = temp1(x(:,ib));
    end
end
   
for iz = 1:size(newRSAmeanH,1)
    delSess(iz,1) = any(isnan(newRSAmeanH(iz,:))); % currently 5 sessions
end

newRSAmeanH = newRSAmeanH(~delSess,:);

% reformat cortex
clear x
% create an index for 
until = size(RSAmean_c_all,1);
x(:,1) = (1:6:until)'; % hit WT
x(:,2) = (4:6:until)';
x(:,3) = (2:6:until)';
x(:,4) = (5:6:until)';
x(:,5) = (3:6:until)';
x(:,6) = (6:6:until)';

newRSAmeanH = [];
for ix = 1:8
    for ib = 1:6
    temp1 = RSAmean_h_all(:,ix);
    newRSAmeanH(:,size(newRSAmeanH,2)+1) = temp1(x(:,ib));
    end
end
   
for iz = 1:size(newRSAmeanH,1)
    delSess(iz,1) = any(isnan(newRSAmeanH(iz,:))); % currently 5 sessions
end

newRSAmeanH = newRSAmeanH(~delSess,:);