mOrange = [[0.901960784313726,0.380392156862745,0.00392156862745098]];
mBlue = [[0.368627450980392,0.235294117647059,0.600000000000000]];


figure(1); clf;
trls = 10;
enc = rand(trls,1);
enc(4) = enc(5) + 2;
enc(10) = enc(10) + 1;
ret = rand(trls,1);
ret(4) = ret(5) + 2;
ret(10) = ret(10) + 0.5;


hold on;
grid on;

b1 = bar([enc, ret]); % Values in the 0-1 range for FaceColor
b1(1).FaceAlpha = 1;
b1(1).FaceColor = mBlue;
b1(2).FaceAlpha = 1;
b1(2).FaceColor = mOrange;
xlabel('Trial Number')
ylabel('Standardized firing rate')
hand = gca;
hand.FontSize = 20;
hand.FontWeight = 'bold';

ewp = enc.*ret;
bEWP = bar(ewp);
bEWP.FaceColor = mean([mBlue;mOrange]);
bEWP.FaceAlpha = 0.2;
ylim([0 10])
xticks([2:2:10])

plot(get(gca,'xlim'), [4 4],'r--','linew',2)

legHand = legend('Encoding', 'Retrieval', 'Encoding x Retrieval');
legHand.FontSize = 16;
legHand.FontWeight = 'normal';
