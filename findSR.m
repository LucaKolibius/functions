function [SR] = findSR
FieldSelection(1) = 1;%timestamps
FieldSelection(2) = 0;
FieldSelection(3) = 0;%sample freq
FieldSelection(4) = 0;
FieldSelection(5) = 1;%samples
ExtractHeader = 1;
ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.

[p2d] = [cd, filesep];
cd('X:\Luca\toolboxes\Neuralynx_19012019\') % after restarting the neuralynx function does not find the path to a mex file, unless you cd into that folder

[CSCfiles] = dir([p2d,'*.ncs']);

dum = cell(1,length(CSCfiles));
it = 1;

[timestampsCSC, dataSamplesCSC,hdrCSC] = Nlx2MatCSC([p2d,CSCfiles(it).name], FieldSelection, ExtractHeader, ExtractMode, []);% import the raw data

% extract the sampling-frequency
chck = regexp(hdrCSC,'SamplingFrequency');
selIdx = [];
for jt = 1:length(chck); selIdx(jt) = ~isempty(chck{jt});end;
selIdx = find(selIdx~=0);
SR = str2double(hdrCSC{selIdx}(min(regexp(hdrCSC{selIdx},'\d')):end));

cd(p2d);
end