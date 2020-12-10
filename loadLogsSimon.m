function [encTrl, retTrl, rtENC, rtRET, missIdx, hitIdx] = loadLogsSimon(path)

%% LOAD IN DATA FROM LOGFILE
opts = delimitedTextImportOptions("NumVariables", 15);

% Specify range and delimiter
opts.DataLines = [1, inf];
opts.Delimiter = ["\t", "_"];

% Specify column names and types
opts.VariableNames = ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"];
opts.SelectedVariableNames = "VarName1";
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "EmptyFieldRule", "auto");


% Import the data
temp = readmatrix(path, opts);

%% Clear temporary variables
clear opts

%% GENERATE ENCTRIAL AND RETTRL
lgfl = str2double(temp);

if ~isnan(lgfl(1)); lgfl = [NaN; lgfl]; end

encTrl = [];
retTrl = [];
encFlag = 0; % ret
counter = 1;
for ii = 2:length(lgfl)
     if isnan(lgfl(ii))
        continue
     end
    
    if isnan(lgfl(ii-1)) && ~isnan(lgfl(ii)) % previous trial nan, this one number
        
        switch encFlag % changes from encoding to retrieval and v.v.
            case 0
                encFlag = 1;
            case 1
                encFlag = 0;
        end
        
    end
    
    switch encFlag
        case 0
            retTrl = [retTrl, counter]; 
        case 1
            encTrl = [encTrl, counter];
    end
    
    counter = counter + 1;
end



%% LOAD IN DATA FROM LOGFILE
opts = delimitedTextImportOptions("NumVariables", 15);

% Specify range and delimiter
opts.DataLines = [1, inf];
opts.Delimiter = ["\t"];

% Specify column names and types
opts.VariableNames = ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"];
opts.SelectedVariableNames = ["VarName1", "Var9"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "EmptyFieldRule", "auto");

% Import the data
temp = readmatrix(path, opts);

isEnc = [];
for ii = 1:size(temp,1)
    
    if strcmp(temp(ii,1),'ENC0')
        saveFlag = 1;
        continue
    elseif strcmp(temp(ii,1), 'RET0')
        saveFlag = 0;
    elseif ~isnan(str2double(temp(ii,1))) % if it is a number
        saveFlag = saveFlag;
    else
        saveFlag = 0;
    end
    
     switch saveFlag
        case 0
            isEnc(ii,1) = 0;
        case 1
            isEnc(ii,1) = 1;
    end
   
end
rtENC = str2double(temp(:,2));
rtENC(~isEnc) = [];




%% LOAD IN DATA FROM LOGFILE
opts = delimitedTextImportOptions("NumVariables", 15);

% Specify range and delimiter
opts.DataLines = [1, inf];
opts.Delimiter = ["\t"];

% Specify column names and types
opts.VariableNames = ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"];
opts.SelectedVariableNames = "Var7";
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15"], "EmptyFieldRule", "auto");

% Import the data
temp = readmatrix(path, opts);

rtRET = str2double(temp);
rtRET(isnan(rtRET)) = [];
rtRET = rtRET(retTrl);

rtRET = [];


% %% RET TRL ORDER
% %% Setup the Import Options and import the data
% opts = delimitedTextImportOptions("NumVariables", 12);
% 
% % Specify range and delimiter
% opts.DataLines = [1, Inf];
% opts.Delimiter = "\t";
% 
% % Specify column names and types
% opts.VariableNames = ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12"];
% opts.SelectedVariableNames = "VarName1";
% opts.VariableTypes = ["double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Specify variable properties
% opts = setvaropts(opts, ["Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12"], "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12"], "EmptyFieldRule", "auto");
% 
% % Import the data
% temp = readtable(path, opts);
% temp =  table2array(temp);
% 
% temp(isnan(temp)) = [];
% retTrlOrd = temp(retTrl);

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 12);

% Specify range and delimiter
opts.DataLines = [1, inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "VarName5", "VarName6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12"];
opts.SelectedVariableNames = ["VarName5", "VarName6"];
opts.VariableTypes = ["string", "string", "string", "string", "double", "double", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12"], "EmptyFieldRule", "auto");

% Import the data
temp = readtable(path, opts);
temp =  table2array(temp);

idx = temp == 1 | temp == 0;
idx = all(idx,2);
temp = temp(idx,:);
miss = sum(temp,2) ~= 2;
hits = sum(temp,2) == 2;

missIdx = find(miss);
hitIdx = find(hits);
end % END OF FUNCTION