% Get input from rpplAnal_pipe and put a breakpointafter calcRppl (around l.76)

filtEEG = microLFP;
for trl = 1:size(microLFP.trial,2)
    for chan = 1:8
        filtEEG.trial{trl}(chan,:) = ft_preproc_bandpassfilter(filtEEG.trial{trl}(chan,:), 1000, [80 140], 3*fix(1000/80)+1, 'fir', 'twopass');
    end
end

counter = 0;
for trl = 1:size(microLFP.trial,2)
    for rip = 1:size(staEnd{trl,1},2)
        counter = counter + 1;
        
        if counter ~= 32
            continue
        else
            disp(trl)
            disp(rip)
        end
        
        %%
        figure(1); clf; hold on;
        mStart = staEnd{trl,1}(rip)-50; if mStart<1; mStart = 1; end;
        mEnd   = staEnd{trl,2}(rip)+50; if mEnd>size(microLFP.trial{trl},2); mEnd = size(microLFP.trial{trl},2); end
        ylabel('\muV')
        
        subplot(211); hold on
        eeg = microLFP.trial{trl}(favChan,mStart:mEnd);
        eeg(1:50) = 5;
        eeg(end-49:end) = 5;
        eeg(104) = 5;
        
        plot(eeg);
        ylim([-5 15])
        xlim([40 length(eeg)-40])
        xticks('')
        ylabel('\muV')
        plot([50 50], get(gca, 'YLim'), 'r--', 'linew', 2)
        plot([length(eeg)-50 length(eeg)-50], get(gca, 'YLim'), 'r--', 'linew', 2)
        %         xticks([50 length(eeg)-50]);
        title('Raw EEG Activity')
        handl = gca;
        handl.FontWeight = 'bold';
        handl.FontSize = 16;
        %%
        subplot(212); hold on;
        eegF = filtEEG.trial{trl}(favChan,mStart:mEnd);
        eegF(1:50) = 0;
        eegF(end-49:end) = 0;
        eegF(104) = 0;
        ylabel('\muV')
        xlabel('Time (ms)')
        plot(eegF);
        %         sgtitle(['Trial: ', num2str(trl), 'ripple# ', num2str(rip)]);
        xlim([40 length(eegF)-40])
        title('EEG activity between 80 and 140 Hz')
        handl = gca;
        handl.FontWeight = 'bold';
        handl.FontSize = 16;
        plot([50 50], get(gca, 'YLim'), 'r--', 'linew', 2)
        plot([length(eeg)-50 length(eeg)-50], get(gca, 'YLim'), 'r--', 'linew', 2)
        
        
        saveas(gcf, ['\\analyse4.psy.gla.ac.uk\project0309\Luca\SFN2020\ripple_example\ripplExample', '_counter',num2str(counter)], 'png'); % emf
        
    end
end