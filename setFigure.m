function setFigure(number)
%% SETUP FIGURE
switch exist('number', 'var')
    case 0
        mFigH = figure('units', 'normalized'); clf;
        
    case 1
        mFigH = figure(number); clf;
        set(mFigH, 'units', 'normalized');
        
end

%     set(gcf, 'units','normalized','outerposition',[0 0 1 1])

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

end % END OF FUNCTION