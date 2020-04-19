function it = manuSorting(dts, it)

cont = input(sprintf('Do you want to continue with wire %d? ', it), 's');
if cont == 'y'
    disp(sprintf('Continuing with wire %d', it));
elseif cont =='n'
    disp('Starting at the first wire');
    it = 1;
end

pnClus = dir('times_*.mat'); % list of all electrode names
for it = it:size(pnClus,1) % loops over all electrode names
%     save('it.mat', 'it');
    disp(['Loading ', pnClus(it).name, ' ', num2str(it), ' / ' , num2str(size(pnClus,1))]);
    figure(1)
    mhandle = wave_clus(pnClus(it).name); % loads results from automatic clustering
    set(mhandle,'WindowStyle','normal'); % undock
%     set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
     set(mhandle,'units','normalized','OuterPosition',[0.9958   -0.0815    1.0083    1.0889]) % fullscreen 
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
    
    % aktivität über die zeit (als eigene funktion)
    disp('Loading longitudinal activity...');
    close all
    mkSpiketimes(dts ,it);
    
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

