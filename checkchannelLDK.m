function checkchannelLDK(ch)
% small helper script that calls Maik's checkunit script on
% all classes in a given times_CSC%i.mat file
% JN 2012-11-26
temp = dir('times_*.mat');
channelsClust={};
for i=1:size(temp,1)
    channelsClust(i) = cellstr(temp(i).name);
end

temp = dir('*_spikes.mat');
channelsDet={};
for i=1:size(temp,1)
    channelsDet(i) = cellstr(temp(i).name);
end

for ch=ch
    channelClust = char(channelsClust(ch));
    channelDet = char(channelsDet(ch));
    
    fname = channelClust(regexp(channelClust,'CSC_')+4:end-4);
    
    if exist(channelClust, 'file')
        load(channelClust);
    else
        warning(['File ' channelClust ' not found!'])
        continue;
    end
    load(channelDet, 'threshold_all');
    
    nc = max(cluster_class(:, 1));
    
    % ?? % delete(sprintf('wav_fig2print_ch%i_cl*.png', channel));
    
    if ~exist('threshold_all','var')
        warning('keine thresholds gespeichert!');
        threshold_all=0;
    end
    
    for i = 1:nc
        if length(find(cluster_class(:,1)==i))>3
            c = checkunit(spikes, ...
                'timevec', (cluster_class(:, 2)-cluster_class(1,2))/1000, ...
                'markvec', cluster_class(:, 1), ...
                'target', i, ...
                'thr', threshold_all);
            
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 1,['\bf wav\_fig2print\_ch' num2str(i) '\_cl' num2str(i)],'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize', 12);
            
            print(c.f, '-r150', ...
                '-dpng', ['wav_fig2print_',fname ,sprintf('_cl%i', i)]);
            close(c.f)
        else
            warning('zu wenig spikes im cluster!');
        end
    end
end


%     % JN 2013-02-25: added windows support
%     if isunix
%         [dummy, hostname] = system('hostname');
%         if isempty(strfind(hostname, 'aachen')) && ...
%                 isempty(strfind(hostname, 'bonn'))
%            system(sprintf('python /media/Projects/johannes/nlx_tools/plot_check.py %i', channel));
%         end
%     elseif ispc
%         system(sprintf('python C:\\Python27\\Lib\\site-packages\\nlx_tools\\plot_check.py %i', channel));
%     end

% cluster correlation. this is based on the function clustercorr from
% treber
%     clustercorr(ch, 100, 200);

channel = ch;
limit = 100;
bins = 200;
maxspikes=5000;
fontSize = 5;
lineWidth = 0.5;

%     load(sprintf('times_CSC%i.mat', channel), 'cluster_class');
numcluster=max(cluster_class(:,1))+1;

for i=1:numcluster
    cluster(i).index_ts=cluster_class(cluster_class(:,1)==i-1,2);
end
h1=figure(800+channel);
set(h1, 'visible','off');
set(h1, 'PaperUnits', 'normalized');
%set(h1, 'PaperSize', 'a4');
set(h1, 'PaperOrientation', 'Portrait')

% this doesn't work anymore under 2018b
%set(gcf, 'PaperUnits', 'Normalized')
%set(gcf, 'PaperSize', [1.2 1.2+max([(numcluster-4) 0])*0.1])
%set(gcf, 'PaperPosition', [0 0 1.2 1.2+max([(numcluster-4) 0])*0.1])
%set(gcf, 'PaperOrientation', 'Portrait')


for i=1:numcluster
    for j=i:numcluster
        subplot(numcluster, numcluster, (i-1)*numcluster+j);
        
        [y maxvalue]=mycorr(cluster(i).index_ts, cluster(j).index_ts, limit/bins, limit);
        y=y./maxvalue;
        bar(-limit:(limit/(bins)):limit, y, 1, 'hist');
        
        title([ 'Cl. ' num2str(i-1) ',' num2str(j-1)], 'FontSize', fontSize);
        %             if (i==1)
        %                 if (i==j)
        %                     title({'Garbage AutoCorr. ', ['Cl. ' num2str(i-1)]}, 'FontSize', fontSize);
        %                 else
        %                     title({'Garbage CrossCorr. ', ['Cl. ' num2str(i-1) ',' num2str(j-1)]}, 'FontSize', fontSize);
        %                 end
        %             elseif (i==j)
        %                 title({'AutoCorr. ', ['Cl. ' num2str(i-1)]}, 'FontSize', fontSize);
        %             else
        %                 title({'CrossCorr. ',[ 'Cl. ' num2str(i-1) ',' num2str(j-1)]}, 'FontSize', fontSize);
        %             end
        xlim([-limit limit]);
        set(gca,'FontSize',fontSize);
        set(gca,'box','on', 'linewidth', lineWidth);
        xlabel('ms', 'FontSize', fontSize);
        ylabel('corr', 'FontSize', fontSize);
    end
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,['\bf Correlations CSC' '\_ch' num2str(channel)],'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize', fontSize*1.2);
%    print(h1,'-dpng', ['ch' num2str(channel) '_Correlations CSC_']);

% overwrite hard-coded ffont size adjustments
%set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal');
%set(findall(0,'type','axes'),'FontSize',fontSize,'fontWeight','normal');
print(h1, '-r150', '-dpng', sprintf('wav_fig2print_ch%i_corr', channel));
close(h1);

end