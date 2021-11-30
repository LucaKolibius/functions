function sanePlot(data, pt)

mFigH = figure(99); clf;
MP = get(0, 'MonitorPositions');
N = size(MP, 1);
% Might want to set an initial position this to some reasonable location
% in the event of the window being "Restored Down".
newPosition = MP(1,:);

if size(MP, 1) == 1
    % Single monitor
    set(mFigH,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
else
    % Multiple monitors - shift to the Nth monitor.
    newPosition(1) = newPosition(1) + MP(N,1);
end
mFigH.set('Position', newPosition, 'units', 'normalized');
mFigH.WindowState = 'maximized'; % Maximize with respect to current monitor.

yfrom = min(data,[], 'all');
yto   = max(data,[], 'all');
for rep = 1:size(data,1)
    plot(data(rep,:), 'color', [0, 0, 0], 'linew', 2); xlim([0 size(data,2)]); ylim([yfrom yto]);
    pause(pt)
end


end % end of function