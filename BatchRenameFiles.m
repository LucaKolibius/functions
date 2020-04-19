% Get all PDF files in the current folder
files = dir('*.ncs');
prev_names=cell(length(files),1);
% Loop through each
for id = 1:length(files)
    prev_names(id) = cellstr(files(id).name);
    % Get the file name (minus the extension)
    [~, f] = fileparts(files(id).name);
    % If numeric, rename
    movefile(files(id).name, sprintf('CSC%d.ncs', id));
end
save('prev_names.mat', 'prev_names');
