% This function loads the SME for all subjects (allSubj)
function [gausConvENC_hitsALL, gausConvENC_missALL, gausConvRET_hitsALL, gausConvRET_missALL] = extractSME(allSubj)
gausConvENC_hitsALL = [];
gausConvENC_missALL = [];
gausConvRET_hitsALL = [];
gausConvRET_missALL = [];
for i = 1 : size(allSubj,1)
subjID = char(allSubj(i));
mSubject = subjID(1:end-3);
mSession = subjID(end-1:end);

try
    cd Z:/Luca/data
catch
    cd /media/ldk898/rds-share/Luca/data
end

if regexp(mSubject, '_')
    mSubject(end) = [];
    mSession = ['S', mSession];
end

cd(mSubject)
cd(mSession)
abc = dir;
cd(abc(3).name)
cd advancedAnalysis
cd SME
SMEdir = cd;
loadVar = dir('SME_*');
load(loadVar.name);

if size(gausConvENC_hits,2)> 1; gausConvENC_hitsALL = [gausConvENC_hitsALL; gausConvENC_hits]; end
if size(gausConvENC_miss,2)>1; gausConvENC_missALL = [gausConvENC_missALL; gausConvENC_miss]; end
% if size(gausConvRET_hits,2)>1; gausConvRET_hitsALL = [gausConvRET_hitsALL; gausConvRET_hits]; end
% if size(gausConvRET_miss,2)>1; gausConvRET_missALL = [gausConvRET_missALL; gausConvRET_miss]; end

end
end
