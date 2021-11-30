%% simulate a bunch of time cells with different list lengths

Listleng=[20 40 60];
t{1,1}=0:0.001:Listleng(1);
t{1,2}=0:0.001:Listleng(2);
t{1,3}=0:0.001:Listleng(3);
timfld=0.3;
nbins=40;
figure;
for n=1:numel(Listleng)
    tfc=round(Listleng(n)*0.3)*1000;
    tfw=2000;
    fr=([zeros(1,tfc-tfw/2) gausswin(tfw)' zeros(1,numel(t{1,n})-tfc-tfw/2)]+0).*10;
    %conv(bb(n,:)',gausswin(200)','same');
    tc{1,n}=P_Neuron(fr,t{1,n});
    %subplot(3,1,n);
    %plot(t{1,n},tc{1,n});
    frconv=conv(tc{1,n},gausswin(200)','same');
    frall(n,:)=[frconv nan(1,(Listleng(3)-Listleng(n))*1000).*3];
    frall_adj(n,:)=imresize(frconv,[1 Listleng(2)*1000]);
    
end
cm=colormap(hot);
figure;imagesc(t{1,3},[1:3],frall,[-0.2 2]);colorbar;colormap(cm);
title('time cells have all different firing fields although they fire proportionally at the same time within each list')

figure;imagesc(t{1,2},[1:3],frall_adj,[-0.2 2]);colorbar;colormap(cm);
title('time cells now fire at the same time within each list')

%% Illustrate c

tuningCurveShuffle = circshift(frconv,randi([1,length(frconv)],1));
figure;plot(t{1,n},frconv, 'Color', [0 0 1]);hold on;
plot(t{1,n},tuningCurveShuffle, 'Color', [1 0 0]);

%% Simulate Poisson Distributed Spikes at certain time point

function [spikeMat]=P_Neuron(fr,tvec)

nbins=length(tvec);
dt=tvec(2)-tvec(1);
for t=1:nbins
    spikeMat(:,t)=double(rand(1,1)<fr(t)*dt);
end

end

%%


