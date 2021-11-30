clear
tic
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\data\allSbj\allPhsPowDat.mat', 'allPhsDat');

allPPC.encNdx = [];
allPPC.encIdx = [];
allPPC.retNdx = [];
allPPC.retIdx = [];
                
for cat = 1:4
    
    for su = 1:length(allPhsDat)
        
        updater = fprintf(1,'Category: %d/%d | SU: %d/%d | timer: %d \n', cat, 4, su, size(allPhsDat,2), toc);

        switch cat
            case 1
                dat = allPhsDat(su).encNdx;
            case 2
                dat = allPhsDat(su).encIdx;
            case 3
                dat = allPhsDat(su).retNdx;
            case 4
                dat = allPhsDat(su).retIdx;
        end
        
        if size(dat,2) < 20
            continue
        end
        
        PPC = [];
        dum = [];
        for curFreq = 1:size(dat,1)
            
            curPhs   = dat(curFreq,:);
            spikeNum = size(dat,2);
            
            for curSpk = 1:size(dat,2)
                dum(curSpk) = nansum( cos ( curPhs(curSpk) - curPhs( curSpk+1:end) ));
            end
            
        PPC(curFreq,1) = nansum(dum) / (spikeNum*(spikeNum-1)/2);
        
        end
        
        switch cat
            
            case 1
                allPPC.encNdx = [allPPC.encNdx PPC];
            case 2
                allPPC.encIdx = [allPPC.encIdx PPC];
            case 3
                allPPC.retNdx = [allPPC.retNdx PPC];
            case 4
                allPPC.retIdx = [allPPC.retIdx PPC];
                
        end
        
        fprintf(repmat('\b', 1, updater));
        
    end % END OF SU LOOP
end % END OF CAT LOOP

figure(1); clf;
ft = 60:2:160;
subplot(221);
plot(ft, mean(allPPC.encIdx,2))
ylim([0 0.3])
title('ENCODING IDX');

subplot(222);
plot(ft, mean(allPPC.encNdx,2))
ylim([0 0.3])
title('ENCODING NDX');

subplot(223);
plot(ft, mean(allPPC.retIdx,2))
ylim([0 0.3])
title('RETRIEVAL IDX');

subplot(224);
plot(ft, mean(allPPC.retNdx,2))
ylim([0 0.3])
title('RETRIEVAL NDX');

sgtitle('PPC averaged over SU with >20 spks in TOI')

%%
dat = [allPhsDat.encNdx];
histogram(dat(10,:))