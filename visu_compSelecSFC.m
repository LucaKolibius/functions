close all
clear

files = dir('X:\Luca\data\engram allocation lfp\segm\sub-*');
filesTrl = dir('X:\Luca\data\engram allocation lfp\trlsegm\sub-*');

nchans = 8;
npad   = 80;
npad2  = npad/2;
srate  = 300;
tv = 1000*(-npad2:npad2)/srate;
getcomp = [0 1 2 3]; % look at the first 4 components
 for ii = 1:2 %: size(files,1)
    load([files(ii).folder, filesep, files(ii).name], 'jdmaps');
    segm = jdmaps;
%     
%     load([filesTrl(ii).folder, filesep, filesTrl(ii).name], 'jdmaps');
%     trlsegm = jdmaps;
    
    pointer = files(ii).name(1:end-4);
    pointer = regexprep(pointer, '_', '-');

    for comp = 1:4
        
        % Continuous Segments
        % component contour
        figure
        subplot(311)
        rmap = reshape(segm(:,end-getcomp(comp))',nchans,npad);
        contourf(npad2+tv(1:end-1),1:nchans,rmap,40,'linecolor','none')
        xlabel(pointer), ylabel('MW Number')
        title(sprintf('Continuous Segments Comp#: %d', comp))

        % TS component
        for cc = 1:8
            subplot(3,8,8+cc)
                        hold on
            plot(npad2+tv(1:end-1),cc+1*zscore(rmap(cc,:)),'k','linew',2)
            plot([40 40], get(gca, 'Ylim'), 'color', 'k');
            xlabel(sprintf('Electrode #%d', cc));
        end
        
        % Powerspectrum
        for cc = 1:8
            subplot(3,8,16+cc)
            plot(linspace(0,srate,200),abs(fft(zscore(rmap(cc,:))/npad,200)).^2,'ks-','linew',2,'markersize',8,'markerfacecolor','w'), hold on
            set(gca,'xlim',[0 40])
            xticks([0:20:100])
        end
        
        
        
%         %% Trial Segments
%         figure
%         subplot(311)
%         title(sprintf('Trial Segments Comp#: %d', comp))
%         rmap = reshape(trlsegm(:,end-getcomp(comp))',nchans,npad);
%         contourf(npad2+tv(1:end-1),1:nchans,rmap,40,'linecolor','none')
%         xlabel(pointer), ylabel('MW Number')
%         title(sprintf('Trial Segments Comp#: %d', comp))
% 
%         
%         % TS component
%         for cc = 1:8
%             subplot(3,8,8+cc)
%             hold on
%             plot(npad2+tv(1:end-1),cc+1*zscore(rmap(cc,:)),'k','linew',2)
%             plot([26 26], get(gca, 'Ylim'), 'color', 'k');
%             xlabel(sprintf('Electrode #%d', cc));
%         end
%         
%         % Powerspectrum
%         for cc = 1:8
%             subplot(3,8,16+cc)
%             plot(linspace(0,srate,200),abs(fft(zscore(rmap(cc,:))/npad,200)).^2,'ks-','linew',2,'markersize',8,'markerfacecolor','w'), hold on
%             set(gca,'xlim',[0 100])
%             xticks([0:20:100])
%         end
    end
end