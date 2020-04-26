%% This script transforms FieldTrip EEG data from clock time into brain time and runs MVPA classification on it.
% By Sander van Bree, María Melcón, and Simon Hanslmayr (April 2020)
% Step 1: Import EEG data, layout and set parameters
% Step 2: Filter data based on parameters fft_min and fft_max
% Step 3: Run ICA on filtered data
% Step 4: Run FFT on components
% Step 5: Find peaks in and phase of components
% Step 6: Continue with components that have peak within frequency range of interest (parameter foi_min and foi_max)
% Step 7: Sort components with most power at peak (descending)
% Step 8: Plot components and let user choose component
% Step 9: Warp component's data to template phase vector (based on peak oscillation)
% Step 10: Run MVPA light on warped data to get brain time TGMs

uiwait(msgbox({'Transforms EEG clock time into EEG brain time and calculates a time generalization matrix';' ';...
    'By Sander van Bree, María Melcón, and Simon Hanslmayr (April 2020)'}))   
uiwait(msgbox({'Required toolboxes: MVPA light, FieldTrip, PsychToolbox';...
    'Required scripts: findpeaks.m';' ';...
    'Before running this script, please ensure that your EEG dataset:';...
    '(1) is in a FieldTrip structure.';...
    '(2) has sufficient data before time of interest (to enable wavelet analysis at early timepoints).';...
    '(3) contains your class 1 and class 2 data, concatenated after each other.';...
    '(4) ideally has the same number of trials for both classes.';...
    '(5) ideally has an even number of trials.'}))

%% Step 1: Import EEG data and layout
%Step 1.1: Import data   
inputpath = inputdlg('Good day. What is the PATH of your FieldTrip formatted EEG dataset? (TO DEV: WRITE maria OR sander');
if strcmp(inputpath,'maria')
    cd('D:\Birmingham\datos\s4')
    load datos_limpios
    cfg         = [];
    cfg.channel = {'all' '-REF' '-HEOG1' '-HEOG2' '-VEOG1' '-VEOG2'};
    data        = ft_preprocessing(cfg,data);
    cd('D:\Birmingham\script')
    layoutfile='biosemi128_1005xx.lay';
elseif strcmp(inputpath,'sander')
    subjnum = inputdlg('Subject name');
    cd \\its-rds\nr-castles\projects\w\wimberm-ieeg-compute\Sander\TGM\ClockToBrainTime\Hippocampal_ClockToBrainTime\SourceData\
    currsubj = strcat('hcdata_',subjnum{1},'.mat');
    load(currsubj)
    data = hcdata;
    cd \\its-rds\nr-castles\projects\w\wimberm-ieeg-compute\Sander\TGM\ClockToBrainTime\
else
    promp     = {'Good day. What is the PATH of your FieldTrip formatted EEG dataset?','What is the NAME of your EEG dataset?','What is the FOLDER of your EEG layout file','What is the NAME of your EEG layout file'};
    answer    = inputdlg(promp);
    inputpath = answer{1};
    inputname = answer{2};
    laypath   = answer{3};
    layname   = answer{4};
    cd(inputpath{1})
    data = load(inputname{1});
    data = data{1};
    cd(laypath{1});
    layoutfile = load(layname{1});
end

%Step 1.2: Parameters
prompt = {'How many components do you want to get from the ICA? (by default, number of sensors)'};
dlgtitle = 'Number of component';
dims = 1;
definput = {num2str(size(data.label,1))};
nc = inputdlg(prompt,dlgtitle,dims,definput);
ncmp  = str2double(nc{1});

promp  = {'To find peaks in your ICA transformed data, we need to apply a band pass filter. What frequency would you like to use for the high pass filter? (e.g. 2)','What frequency would you like to use for the low pass filter? (e.g. 30)'};
fft  = inputdlg(promp);
fft_min  = str2double(fft{1});
fft_max  = str2double(fft{2});

promp  = {'We later filter ICA components to include only components with peaks at your frequency range of interest. What is the lower bound frequency of interest? (e.g. 8)','What is the upper bound frequency of interest? (e.g. 12)'};
freq  = inputdlg(promp);
freq_min  = str2double(freq{1});
freq_max  = str2double(freq{2});

% promp  = {'Resampling will speed up the script a lot. What should be the new sampling rate? (e.g. 256; cancel = no resampling)'};
% resample  = inputdlg(promp);
% if isempty(resample{1})~=1
% cfg = [];
% cfg.resamplefs = str2double(resample{1});
% data = ft_resampledata(cfg,data);
% end

promp  = {'To perform the time frequency analysis, what is your time window of interest on which it should be centered (in second)? We use 0.005 s steps. From (e.g. 0)','To (e.g. 1)'};
toi  = inputdlg(promp);
t_min  = str2double(toi{1});
t_max  = str2double(toi{2});

toi_min  = data.time{1}(1,findnearest(data.time{1},t_min));
toi_max  = data.time{1}(1,findnearest(data.time{1},t_max));
toi_step = data.time{1}(1,2)-data.time{1}(1,1);% original rata sampling rate for frequency analysis

%% Step 2: Filter data
cfg               = [];
cfg.bpfilter      = 'yes';
cfg.bpfreq        = [fft_min fft_max]; %based on input parameters
data_bp           = ft_preprocessing(cfg,data);

%% Step 3: ICA
cfg               = [];
cfg.method        = 'runica';
cfg.runica.pca    = ncmp; %obtain N component, to reduce time
comp              = ft_componentanalysis(cfg ,data_bp);

%Create template data with full parameter info to run FFT on
temp_comp=data;
temp_comp.trial=comp.trial;
temp_comp.label={};
for cmp=1:length(comp.label)
    temp_comp.label{cmp,1}=data.label{cmp,1};
end

%% Step 4: Perform frequency analysis on the components and calculate phase
cfgtf           = [];
cfgtf.method    = 'wavelet';
cfgtf.width     = 5;
cfgtf.toi       = toi_min:toi_step:toi_max;
cfgtf.foi       = (fft_min:1:fft_max);
cfgtf.output    = 'fourier';
fspec           = ft_freqanalysis(cfgtf,temp_comp);
fspec.powspctrm = abs(fspec.fourierspctrm);

%% Step 5: Find peaks in and phase of components
for cmp = 1:ncmp
    powtf(:,:,cmp)=squeeze(nanmean(fspec.powspctrm(:,cmp,:,:),1));
    pspec(:,cmp)=squeeze(nanmean(powtf(:,:,cmp),2)); %get power spectrum of components
    [~,oscpeak{cmp}]=findpeaks(pspec(:,cmp),'SortStr','descend');% what's the dominant oscillation for current chan?
    if isempty(oscpeak{cmp}) == 1 %if no clear peak, grab max power channel
        oscpeak{cmp} = find(pspec(:,cmp)==max(pspec(:,cmp)));
    else %ordinarily just pick the strongest peak
        oscpeak{cmp} = oscpeak{cmp}(1);
    end
    phs{cmp}=squeeze(angle(fspec.fourierspctrm(:,cmp,oscpeak{cmp},:))); %Simon: don't average around dominant osc, just use dominant osc
end

oscpeak = cell2mat(oscpeak); %get the dominant oscillation vector into array

%% Step 6: Find components that have a dominant oscillation in the frequency of interest range
foicomps = [];
while isempty(foicomps)
    foicomps = find(fspec.freq(oscpeak)>=freq_min & fspec.freq(oscpeak)<=freq_max); %grab only those comps
    if isempty(foicomps)
        foi_range = inputdlg('WARNING: None of the comps shows a peak oscillation at your frequency of interest. Please re-enter a new minimum and maximum frequency of interest separated by a comma (e.g: 7,12)');
        freq_min = (extractBefore(foi_range,","));
        freq_min = str2double(freq_min{1});
        freq_max = (extractAfter(foi_range,","));
        freq_max = str2double(freq_max{1});
    end
end

ind1 = 1;
pspec_comp = zeros(numel(foicomps),2);
for ind2 = foicomps
    pspec_comp(ind1,1)=pspec(oscpeak(ind2),ind2);
    pspec_comp(ind1,2)=ind2;
    ind1 = ind1+1;
end
pspec_comp = round(pspec_comp);

%% Step 7: Sort components with most power at peak (descending)
pspec_comp=sortrows(pspec_comp,1,'descend'); %sort components from highest to lowest average power

%% Step 8: Plot components and let user choose component
figure
coi = [];%component of interest
count=1;
cxs_max=max(max(max(powtf(oscpeak(foicomps),:,foicomps))));
ic=1;

uiwait(msgbox({'Please decide one component, factoring in wheter the component:';' ';...
    '(1) has a topography of interest.';...
    '(2) has high power at the frequency of interest.';...
    '(3) optionally, the variance explained by the component (lower component number means higher r^2).'}))
uiwait(msgbox({'Instructions to browse through components:';' ';...
    'Press forward/back arrow to see the next/previous component';...
    'Once you have decided for one component, click in that component and press Q to quite visualization'}))
while ic <= size(pspec_comp,1)
    % component topography
    subplot(5,2,[1 3 5 7 9]);
    cfg           = [];
    cfg.component = pspec_comp(ic,2); % specify the component(s) that should be plotted
    cfg.layout    = layoutfile; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, comp)
    colorbar
    lim=max(abs(comp.topo(:,pspec_comp(ic,2))));
    lim=lim+(lim/100*10);
    caxis([-lim lim])
    title({[num2str(ic) '/' num2str(size(pspec_comp,1)) ' Components'] ['Component ' num2str(pspec_comp(ic,2))]})

    % time-frequency plot
    subplot(5,2,[2 4 6 8]);
    pcolor(fspec.time,fspec.freq,powtf(:,:,pspec_comp(ic,2)));
    shading interp
    ylim([fft_min fft_max]);
    caxis([0 cxs_max])
    title(sprintf('Power: %s',pspec_comp(ic,1)))
    
    % dominant oscillation plot
    subplot(5,2,10);
    plot(squeeze(pspec(:,pspec_comp(ic,2),1)),fspec.freq);
    ylim([fft_min fft_max]);
    xlim([0 cxs_max])
    
    keydown = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    [keyIsDown,timeStamp,keyCode] = KbCheck;
    key = KbName(find(keyCode));
    if (keydown == 0) %grab the component if click
        coi(count)=pspec_comp(ic,2);
        count=count+1;
    elseif value == 28 %go back previous component if press back arrow
        disp(value);
        ic = ic-2;
    elseif strcmp(key,'q') %stop the loop if it is not necessesary to keep visualising
        close
        break
    end
    ic = ic+1;
end

%remove the component from original data
cfg           = [];
cfg.component = coi;
data          = ft_rejectcomponent (cfg, comp, data_bp);

% cut out the time window of interest
cfg        = [];
cfg.toilim = [t_min t_max];
data       = ft_redefinetrial(cfg, data);

%% Step 9: Warp component's data to template phase vector (based on peak oscillation)
% Re-Organize EEG data by phase
% split data into two datasets by condition
ntrls=length(data.trial);
cfg=[];
cfg.trials=1:round(ntrls/2);
dat1=ft_selectdata(cfg,data);
cfg.trials=round(ntrls/2)+1:ntrls;
dat2=ft_selectdata(cfg,data);

%Get phase for each condition
dat1_bt=dat1;
dat2_bt=dat2;
phs1=phs{coi}(1:round(ntrls/2),:);
phs2=phs{coi}(round(ntrls/2)+1:end,:);
nsec=data.time{1}(1,end)-data.time{1}(1,1);
Ncycles=fspec.freq(oscpeak(coi))*nsec; %number of cycles * seconds
promp  = {'To dynamically time warp the original data based on the chosen component''s data, we need a new sampling rate. High sampling rate is important for high frequencies of interest, but slows down classificaiton a lot. (e.g. 256.)'};
phs_sr = inputdlg(promp);
phs_sr = str2double(phs_sr{1});
tempphs=linspace(-pi,(2*pi*Ncycles)-pi,phs_sr);% set up phase bins for unwrapped phase (angular frequency)
timephs=linspace(0,Ncycles,phs_sr)+1; %time vector of the unwrapper phase

for nt=1:round(ntrls/2)
    tmpphs1=unwrap(phs1(nt,:));
    tmpphs2=unwrap(phs2(nt,:));
    % Warp phase of single trial onto template phase
    [~,ix1,~] = dtw(tmpphs1,tempphs);
    [~,ix2,~] = dtw(tmpphs2,tempphs);

    % Create warped trials
    tmptrl1=dat1.trial{1,nt}(:,ix1);
    tmptrl2=dat2.trial{1,nt}(:,ix2);
    dat1_bt.trial{1,nt}=imresize(tmptrl1,[size(tmptrl1,1) numel(tempphs)]);
    dat2_bt.trial{1,nt}=imresize(tmptrl2,[size(tmptrl2,1) numel(tempphs)]);    
    dat1_bt.time{1,nt}=timephs;
    dat2_bt.time{1,nt}=timephs;
end


%% Step 10: Run MVPA light on warped data to get brain time TGMs

%append data into one set
cfg      = [];
fulldata = ft_appenddata(cfg,dat1_bt,dat2_bt);

%append data into one set
cfg = [];
cfg.keeptrials = 'yes';
cfg.removemean = 'no';
cfg.vartrllength = 2;
fulldata = ft_timelockanalysis(cfg,fulldata);

%define clabels
clabel1 = ones(round(ntrls/2),1);
clabel2 = 2*ones(round(ntrls/2),1);

% Put them together in one dataset
clabel = [clabel1;clabel2];

% Run MVPA light
cfg =  [];
zscoreopt = inputdlg('Would you like to perform global z-scoring of the data before running MVPA classification? Type y or n');
if strcmp(zscoreopt,'y')
    cfg.preprocess  = 'zscore';
end
perfmetric = inputdlg('Would you like the performance metric for classification to be decision-values (dval; distance from hyperplane) or accuracy (acc)? write acc or dval');
cfg.classifier  = 'lda';
cfg.metric      = perfmetric{1};
cfg.repeat      = 5;
cfg.cv          = 'kfold';
cfg.k           = 5;
[TGM_perf, TGM_perf_specs] = mv_classify_timextime(cfg, fulldata.trial, clabel);

%if user chose dval, continue with only one class
if strcmp(perfmetric{1},'dval')==1
TGM_perf = squeeze(TGM_perf(:,1,:));
end

zscoreopt2 = inputdlg('Would you like to perform z-scoring of the TGMs? Recommended for decision value data. Type y or n');
if strcmp(zscoreopt2,'y')
   TGM_perf = zscore(TGM_perf,0,'all');
end

%Grab time vector
timevec = fulldata.time;

% Plot brain time TGMs
smoothlevel = inputdlg('What should be the sd of the Gaussian smoothing kernel for the TGMs? (default: 2)');
if isempty(smoothlevel{1})
    smoothlevel = 2;
else
    smoothlevel = str2double(smoothlevel{1});
end

figure
hold on
cfg_plot= [];
cfg_plot.x   = timevec;
cfg_plot.y   = cfg_plot.x;
TGM_perf_smooth = imgaussfilt(TGM_perf,smoothlevel);
mv_plot_2D(cfg_plot, TGM_perf_smooth);
xlim([timevec(1) timevec(end)]);
ylim([timevec(1) timevec(end)]);
xticks(yticks) % make ticks the same on the two axes
%   colormap jet
title('Average accuracy sorted by brain time')
ylabel('Cycle number training data')
xlabel('Cycle number test data')

%% Step 11: Save TGM
savefold = strcat(savepath,'TGM_braintime');
save(savefold, 'TGM_perf_smooth', 'TGM_perf_specs', 'TGM_perf', 'clabel', 'timevec', '-v7.3');
