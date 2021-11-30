function PPC_linenoise(data, idxSpks)
npnts   = 3000; % trial length
sr = 1000;
trls    = floor(size(data.lfp_noArt,2)/npnts);
dt = 0:npnts;
nchan = 8;
for chan = 1 : nchan

for trl = 1 : trls
    templfp =  data.lfp_noArt(chan,npnts*(trl-1)+1:npnts*trl);
    tempspks = idxSpks(idxSpks>=npnts*(trl-1)+1 & idxSpks<=npnts*trl) - npnts*(trl-1);
    tempspks = histcounts(tempspks,dt);
    
    data_ppc.sampleinfo(trl,1) = npnts*(trl-1)+1;
    data_ppc.sampleinfo(trl,2) = npnts*trl;
    
    data_ppc.label = {'EEG'; 'Spks'};
    data_ppc.trial{1,trl} = [templfp; tempspks];
    data_ppc.time{1,trl} = 1/sr :1/sr : 3;
end

    cfg=[];
    cfg.channel=1;
    EEGdata=ft_selectdata(cfg,data_ppc);
    
    cfgtf=[];
    cfgtf.method='wavelet';
    cfgtf.width=8;
    cfgtf.toi=0.001:0.001:3;
    cfgtf.foi=[10:2:150];
    cfgtf.output='fourier';
    fspec=ft_freqanalysis(cfgtf,EEGdata);
    fspec.powspctrm=abs(fspec.fourierspctrm);
    powtf=squeeze(nanmean(fspec.powspctrm,1));
    ps=squeeze(nanmean(powtf,2));
    [tmpf,tmpl]=findpeaks(ps,'SortStr','descend');% find dominant oscillation
    
    % Extract complex values for PPC calculation at spike times
    phi=squeeze(fspec.fourierspctrm(:,1,:,:));% extract phase
    ntrls=length(data_ppc.trial);
    tofi(1)=find(data_ppc.time{1,1}==cfgtf.toi(1));
    tofi(2)=find(data_ppc.time{1,1}==cfgtf.toi(end));
    for nt=1:ntrls
        spkts(nt,:)=data_ppc.trial{1,nt}(2,tofi(1):tofi(2));
    end
    
    [ tmp ] = phi;
    [ ts ]  = spkts;
    ix = find(ts);
    tmp1 = permute(phi,[1 3 2]);
    tmp2 = reshape(tmp1,[size(tmp1,1)*size(tmp1,2) size(tmp1,3)]);
    [ PPC ] = computeSPK2LFPcoupling_S( tmp2(ix,:), 'ppc');
    
    
    %% Plot Powerspectrum of simulated Signal, just to check it looks EEG'ish;
    % You can comment this whole section out once you are satisfied with the
    % EEG signal
    figure;
    subplot(1,10,1:6);
    pcolor(fspec.time,fspec.freq,squeeze(powtf));
    shading interp
    ylim([10 150]);
    hold on
    subplot(1,10,7:8);
    plot(ps,fspec.freq);
    ylim([10 150]);
    hold off
    subplot(1,10,9:10);
    plot(PPC,fspec.freq);
    ylim([10 150]);
    hold off
    
    cd('X:\Luca\visu\ppc_linenoise')
saveas(gcf, sprintf('PPCchan%d_%d', chan, round(toc)), 'png');
end
end
