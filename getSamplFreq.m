function [FS] = getSamplFreq(p2d)

FieldSelection(1) = 1;%timestamps
FieldSelection(2) = 0;
FieldSelection(3) = 0;%sample freq
FieldSelection(4) = 0;
FieldSelection(5) = 1;%samples
ExtractHeader = 1;
ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.

[CSCfiles] = dir([p2d,'*.ncs']);
it = 1;

[~, ~,hdrCSC] = Nlx2MatCSC([CSCfiles(it).folder, filesep, CSCfiles(it).name], FieldSelection, ExtractHeader, ExtractMode, []);% import the raw data

%extract the scale factor
chck   = regexp(hdrCSC,'ADBitVolts');
selIdx = [];
for jt = 1:length(chck); selIdx(jt) = ~isempty(chck{jt});end;
selIdx = find(selIdx~=0);
scalef = str2double(hdrCSC{selIdx}(min(regexp(hdrCSC{selIdx},'\d')):end));

% extract the sampling-frequency
chck = regexp(hdrCSC,'SamplingFrequency');
selIdx = [];
for jt = 1:length(chck); selIdx(jt) = ~isempty(chck{jt});end;
selIdx = find(selIdx~=0);
FS = str2double(hdrCSC{selIdx}(min(regexp(hdrCSC{selIdx},'\d')):end));

end