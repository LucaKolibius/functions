function spkRppl_visu(resQuant, resCoinc, resTime, sets)% general settings
coinc   = sets.coinc;   % 0 = QUANTITY | 1 = COINCIDENCE
avgTrl  = sets.avgTrl;  % 0 for each trial independently | 1 for average over one SU
del0    = sets.del0;    % 0 = keep 0s  | 1 = delete 0s (this isn't implemented for avgSU - either kick out the trials with 0 or the SU with 0...)
medmean = sets.medmean; % 0 for mean | 1 for median
tw      = sets.tw;

% % 
if coinc == 0
    dat = resQuant;
    dt = 0:2:50;
elseif coinc == 1
    dat = resCoinc;
    dt = 0:1:25;
end

% PREPRO DATA
% PER TRIAL
if avgTrl == 0
    allBund = vertcat(dat{:});
    
    if del0 == 1
        allBund(allBund(:,1) == 0,:) = []; % delete 0s
    end
    
    IUidx = allBund(allBund(:,2) == 1, 1);
    IUndx = allBund(allBund(:,2) == 2, 1);
    GUidx = allBund(allBund(:,2) == 3, 1);
    GUndx = allBund(allBund(:,2) == 4, 1);
    SUnbd = allBund(allBund(:,2) == 6, 1);
    
    if coinc == 1
        SUbnd = allBund(allBund(:,2) == 5, 1);
    end
    
% AVERAGED OVER SU
elseif avgTrl == 1
   
    IUidx = [];
    IUndx = [];
    GUidx = [];
    GUndx = [];
    SUnbd = [];
    
    for bund = 1 : size(dat,2)
        if isempty(dat{bund})
            continue
        end
        
        temp = dat{bund}(dat{bund}(:,2) == 1, 1);
        IUidx = [IUidx; mean(temp)];
        
        temp = dat{bund}(dat{bund}(:,2) == 2, 1);
        IUndx = [IUndx; mean(temp)];
        
        temp = dat{bund}(dat{bund}(:,2) == 3, 1);
        GUidx = [GUidx; mean(temp)];
        
        temp = dat{bund}(dat{bund}(:,2) == 4, 1);
        GUndx = [GUndx; mean(temp)];
        
        if coinc == 1
            temp = resCoinc{bund}(resCoinc{bund}(:,2) == 5, 1);
            SUbnd = [SUbnd; nanmean(temp)];
        end
        
        temp = dat{bund}(dat{bund}(:,2) == 6, 1);
        SUnbd = [SUnbd; mean(temp)];
    end
end

%  ACTUAL VISUALISATION
%  INDEX UNITS
close all
figure('units','normalized','outerposition',[0 0 1 1])

% SUPER-TITLE
tits = '';
if coinc == 0
    tits = 'Ripple Quantity |';
elseif coinc == 1
    tits = 'SPK-RPPL Coincidences |';
end

if tw == 0
    tits = strjoin({tits, 'Whole Trial |'});
elseif tw == 1
    tits = strjoin({tits, 'Pre-Cue |'});
elseif tw == 2
    tits = strjoin({tits, 'Peri-Cue |'});
elseif tw == 3
    tits = strjoin({tits, 'Peri-Resp |'});
end

if avgTrl == 1
    tits = strjoin({tits, 'Average per SU'});
else
    tits = strjoin({tits, 'Every Trial'});
end

if del0 == 1
    tits = strjoin({tits, '| No Zeros'});
end

sgtitle(tits)


subplot(311)
hold on

handl = histogram(IUidx, dt, 'Normalization','probability');
mhist = histogram(IUndx, dt, 'Normalization','probability');
mhist.FaceColor = [0.93,0.69,0.13];

switch medmean
    case 0 % mean
        plot([nanmean(IUndx) nanmean(IUndx)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        plot([nanmean(IUndx) nanmean(IUndx)], get(gca, 'YLim'), 'color', [0.93,0.69,0.13], 'linew', 2, 'linestyle', '-');
        plot([nanmean(IUidx) nanmean(IUidx)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        title(sprintf('Index Units - Indexed Trials (#%d; M = %.2f) vs. Non-Indexed Trials (#%d; M = %.2f)', size(IUidx,1)-sum(isnan(IUidx)), nanmean(IUidx), size(IUndx,1)-sum(isnan(IUndx)), nanmean(IUndx)))
    case 1 % median
        plot([nanmedian(IUndx) nanmedian(IUndx)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        plot([nanmedian(IUndx) nanmedian(IUndx)], get(gca, 'YLim'), 'color', [0.93,0.69,0.13], 'linew', 2, 'linestyle', '-');
        plot([nanmedian(IUidx) nanmedian(IUidx)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        title(sprintf('Index Units - Indexed Trials (#%d; MED = %.2f) vs. Non-Indexed Trials (#%d; MED = %.2f)', size(IUidx,1)-sum(isnan(IUidx)),nanmedian(IUidx), size(IUndx,1)-sum(isnan(IUndx)), nanmedian(IUndx)))
end

% if coinc == 0
%     xlabel('Quantitiy of Ripples');
% elseif coinc == 1
%     xlabel('Spike-Ripple Coincidences');
% end

ylabel('Percentage')
legend([handl, mhist], 'Indexed', 'Non-Indexed')

% GRAY UNITS
subplot(312)
hold on
histogram(GUidx, dt, 'Normalization','probability');
mhist = histogram(GUndx, dt, 'Normalization','probability');
mhist.FaceColor = [0.93,0.69,0.13];

switch medmean
    case 0
        plot([nanmean(GUndx) nanmean(GUndx)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        plot([nanmean(GUndx) nanmean(GUndx)], get(gca, 'YLim'), 'color', [0.93,0.69,0.13], 'linew', 3, 'linestyle', '-');
        plot([nanmean(GUidx) nanmean(GUidx)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        title(sprintf('Gray Units - Indexed Trials (#%d; M = %.2f) vs. Non-Indexed Trials (#%d; M = %.2f)', size(GUidx,1)-sum(isnan(GUidx)), nanmean(GUidx), size(GUndx,1)-sum(isnan(GUndx)), nanmean(GUndx)))
    case 1
        plot([nanmedian(GUndx) nanmedian(GUndx)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        plot([nanmedian(GUndx) nanmedian(GUndx)], get(gca, 'YLim'), 'color', [0.93,0.69,0.13], 'linew', 3, 'linestyle', '-');
        plot([nanmedian(GUidx) nanmedian(GUidx)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        title(sprintf('Gray Units - Indexed Trials (#%d; MED = %.2f) vs. Non-Indexed Trials (#%d; MED = %.2f)', size(GUidx,1)-sum(isnan(GUidx)), nanmedian(GUidx), size(GUndx,1)-sum(isnan(GUndx)), nanmedian(GUndx)))
end


% SINGLE UNITS
subplot(313)
hold on
mhist = histogram(SUnbd, dt, 'Normalization','probability');
mhist.FaceColor = [0.93,0.69,0.13];

if coinc == 1
    handl = histogram(SUbnd, dt, 'Normalization','probability');
    handl.FaceColor = [0.00,0.45,0.74];
end

switch medmean
    case 0
        plot([nanmean(SUnbd) nanmean(SUnbd)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        plot([nanmean(SUnbd) nanmean(SUnbd)], get(gca, 'YLim'), 'color', [0.93,0.69,0.13], 'linew', 2.5, 'linestyle', '-');
        title(sprintf('Single Units (M = %.2f))', nanmean(SUnbd)))
        if coinc == 1
            plot([nanmean(SUbnd) nanmean(SUbnd)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
            title(sprintf('Single Units - On-Bundle (M = %.2f) vs. Off-Bundle (M = %.2f)', nanmean(SUbnd), nanmean(SUnbd)))
        end
    case 1
        plot([nanmedian(SUnbd) nanmedian(SUnbd)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        plot([nanmedian(SUnbd) nanmedian(SUnbd)], get(gca, 'YLim'), 'color', [0.93,0.69,0.13], 'linew', 2.5, 'linestyle', '-');
        title(sprintf('Single Units (MED = %.2f))', nanmedian(SUnbd)))
        if coinc == 1
            plot([nanmedian(SUnbd) nanmedian(SUnbd)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
            title(sprintf('Single Units - On-Bundle (MED = %.2f) vs. Off-Bundle (MED = %.2f)', nanmedian(SUbnd), nanmedian(SUnbd)))
        end
end

if coinc == 1
    legend([handl, mhist], 'On-Bundle', 'Off-Bundle')
end

% within SU avg
if avgTrl == 1
try
GUdiff = GUidx - GUndx;
% close all
figure('units','normalized','outerposition',[0 0 1 1])
subplot(10,1,2:10)
sgtitle(tits)
hold on
dt = min(GUdiff):1:max(GUdiff);
handl = histogram(GUdiff, dt, 'Normalization','probability');

switch medmean
    case 0 % mean
        plot([nanmean(GUdiff) nanmean(GUdiff)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        plot([nanmean(GUdiff) nanmean(GUdiff)], get(gca, 'YLim'), 'color', [0.93,0.69,0.13], 'linew', 2, 'linestyle', '-');
        plot([nanmean(GUdiff) nanmean(GUdiff)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        title(sprintf('Gray Units: Idx-Ndx (M = %.2f)', nanmean(GUdiff)))
    case 1 % median
        plot([nanmedian(GUdiff) nanmedian(GUdiff)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        plot([nanmedian(GUdiff) nanmedian(GUdiff)], get(gca, 'YLim'), 'color', [0.93,0.69,0.13], 'linew', 2, 'linestyle', '-');
        plot([nanmedian(GUdiff) nanmedian(GUdiff)], get(gca, 'YLim'), 'color', 'k', 'linew', 4, 'linestyle', '-');
        title(sprintf('Gray Units: Idx-Ndx (MED = %.2f)', nanmedian(GUdiff)))
end


catch
end % end of try catch
end % end of if condition

%% TIMING VISU
% if avgTrl == 1 && tw > 1
if tw > 1
%  PER TRIAL
allBund = vertcat(resTime{:});                           % get all data
% allBund(allBund(:,1) == 0,:) = [];                     % delete 0s

allIndx = vertcat(resQuant{:}); allIndx = allIndx(:,2);  % get identifier (1-6)
IUidx = mean(allBund(allIndx == 1, :),1);
IUndx = mean(allBund(allIndx == 2, :),1);
GUidx = mean(allBund(allIndx == 3, :),1);
GUndx = mean(allBund(allIndx == 4, :),1);
SUnbd = mean(allBund(allIndx == 6, :),1);

dt = -50:100:2001;
dt(1) = 1;
% dt = 1:1:2001;
TP = zeros(6,size(dt,2)-1);
for ii = 1 : size(dt,2)-1

    TP(1,ii) =  sum(IUidx(dt(ii):dt(ii+1)));
    TP(2,ii) =  sum(IUndx(dt(ii):dt(ii+1)));
    TP(3,ii) =  sum(GUidx(dt(ii):dt(ii+1)));
    TP(4,ii) =  sum(GUndx(dt(ii):dt(ii+1)));
    TP(5,ii) =  sum(SUnbd(dt(ii):dt(ii+1)));

end

% maybe delete first and last entry of TP instead of this stupid xlim business

tits2 = {'Index Units - Indexed', 'Index Units - NonIndexed', 'Gray Units - Indexed', 'Gray Units - NonIndexed', 'Single Units'};
figure('units','normalized','outerposition',[0 0 1 1])
subplot(10,1,2:10)


sgtitle(tits)

for pp = 1:5
subplot(5,1,pp)
plot(TP(pp,:))
title(tits2(pp))
ylim([0 3])
xticks([])
hold on
plot([11 11], [0 3], 'linew', 3)
xlim([2 20])
end

subplot(511)
xticks([1:1:20]);
xticklabels([-1000:100:1000])
xlabel('Time to Resp [ms]')
ylabel('Ripple Activity')

end

end % end of function