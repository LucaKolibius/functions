% ideas:
% segmentierung von spike detection sollte auf varianz basieren!

clear all
addpath(genpath('Z:\Common\fieldtrip-20170618')); % fieldtrip
ft_defaults;
addpath(genpath('Z:\Luca\functions\wave_clus-master')); % waveclus 3.0
addpath('Z:\Luca\functions'); % my functions
addpath('Z:\Luca\functions\Neuralynx_19012019'); % Neuralynx (the commons folder function doesnt work on my PC)
addpath(genpath('Z:\Luca\TREBER\Scripts'))

try
p2d = cd;
p2d(end+1)='\';
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, hitsIdx, missIdx, allTrials, sessionNum, retTrigger, encTrigger]=loadLogs(p2d);
dts.start = encTrigger(1);
dts.end = max(max(retTrigger))+1; % remember that retTrigger is not linearly increasing
dts.dt = 0.001:0.001:dt_end; % create a linerally spaced vector from the beginning to the end of the experiment in steps of 1ms

catch
    disp('Using default dt values');
    dts.start = 1;
    dts.end = 3000;
    dts.dt = 3;
end

% positive detection
par = set_parameters();
par.detection = 'pos';
Get_spikes('all', 'parallel', true, 'par', par);
movefile('*.mat', 'posDetect'); % moves all .mat files that include the just detected spikes into the folder "posDetect"

% % negative detection
par.detection = 'neg';
Get_spikes('all', 'parallel', true, 'par', par);
movefile('*.mat', 'negDetect'); % moves all .mat files that include the just detected spikes into the folder "negDetect"

% positive clustering
cd posDetect
Do_clustering('all', 'parallel', true, 'make_plots', true);
posClus = dir('times_*mat');

% % negative clustering
cd ..\negDetect
Do_clustering('all', 'parallel', true, 'make_plots', true);
negClus = dir('times_*mat');

clear it
% positive manual inspection
cd ..\posDetect
if ~exist('it')
    it=1;
end
it = manuSorting(it, dts);

%     checkchannelLDK(it)
%     mkSpiketimes(dts ,it);

% negative manual inspection
cd ..\negDetect
if ~exist('it')
    it=1;
end
it = manuSorting(it, dts);

for i = i:size(negClus,1)
    disp(['loading ', negClus(i).name, ' ', num2str(i), ' / ' , num2str(size(negClus,1))]);
    figure(1)
    mhandle = wave_clus(negClus(i).name); % loads results from automatic clustering
    set(mhandle,'WindowStyle','normal'); % undock
    set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
    
    % wait until you enter 'y'
    nxt = input('Load Checkchannel? ', 's');
    while nxt~='y' && nxt~='n'
        nxt = input('Load Checkchannel? ', 's');
    end
    
    
    if nxt=='y'
        disp('loading checkchannel');
        figure(2);
        checkchannelLDK(i)
    elseif nxt=='n'
        continue
    end
    
    % aktivität über die zeit (als eigene funktion)
    disp('loading longitudinal activity');
    close all
    mkSpiketimes(dt_start,dt_end,dt,i);
    
    nxt = input('Do you want to reopen waveclus? ', 's');
    while nxt~='y' && nxt~='n' % 121 = 'y' // 110 = 'n'
        nxt = input('Do you want to reopen waveclus? ', 's');
    end
    
    if nxt=='y'
        figure(1)
        mhandle = wave_clus(negClus(i).name);
        set(mhandle,'WindowStyle','normal'); % undock
        set(mhandle,'units','normalized','OuterPosition',[-0.004 0.03 1.008 0.978], 'InnerPosition', [0 0.037 1 0.892]) % fullscreen
        
        while nxt~='y' && nxt~='n' % 121 = 'y'
            nxt = input('Do you want to reopen waveclus? ', 's');
        end
    elseif nxt=='n'
        continue
    end
end
clear i

