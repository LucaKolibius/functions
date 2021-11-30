function it = manuSorting(it)

%% SETUP FIGURE
mFigH = figure(1); clf;

MP = get(0, 'MonitorPositions');
N = size(MP, 1);

% Might want to set an initial position this to some reasonable location
% in the event of the window being "Restored Down".
newPosition = MP(1,:);

% Multiple monitors - shift to the Nth monitor.
newPosition(1) = newPosition(1) + MP(N,1);
mFigH.set('Position', newPosition, 'units', 'normalized');
mFigH.WindowState = 'maximized'; % Maximize with respect to current monitor.


cont = input(sprintf('Do you want to continue with wire %d? ', it), 's');
if cont == 'y'
    disp(sprintf('Continuing with wire %d', it));
elseif cont =='n'
    disp('Starting at the first wire');
    it = 1;
end

pnClus = dir('times_*.mat'); % list of all electrode names
for it = it:size(pnClus,1) % loops over all electrode names
    disp(['Loading ', pnClus(it).name, ' ', num2str(it), ' / ' , num2str(size(pnClus,1))]);
    figure(1)
    mhandle = wave_clus(pnClus(it).name); % loads results from automatic clustering
    set(mhandle,'WindowStyle','normal'); % undock
    
    set(mhandle,'units','normalized','OuterPosition', newPosition) % fullscreen
    figure(1)
    set(mhandle, 'WindowState', 'maximized');
    close(1);
     
    % wait until you enter 'y' or 'n'
    nxt = input('Load Checkchannel? ', 's');
    while nxt~='y' && nxt~='n'
        nxt = input('Load Checkchannel? ', 's');
    end
    
    
    if nxt=='y'
        disp('Loading checkchannel...');
        figure(2);
        checkchannelLDK(it)
    elseif nxt=='n'
        continue
    end
    
    % aktivit�t �ber die zeit (als eigene funktion)
    disp('Loading longitudinal activity...');
    close all
%     mkSpiketimes(dts ,it);
    
    nxt = input('Do you want to reopen waveclus? ', 's');
    while nxt~='y' && nxt~='n' % 121 = 'y' // 110 = 'n'
        nxt = input('Do you want to reopen waveclus? ', 's');
    end
    
    if nxt=='y'
        figure(1)
        mhandle = wave_clus(pnClus(it).name);
        set(mhandle,'WindowStyle','normal'); % undock
        set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
        close(1)
        nxt = input('Load next wire? ', 's');
        while nxt~='y' && nxt~='n'
            nxt = input('Load next wire? ', 's');
        end
        if nxt == 'y'
            continue
        elseif nxt == 'n'
            disp('Loading checkchannel...');
            checkchannelLDK(it)
        end
        
    elseif nxt=='n'
        continue
    end
end
end

