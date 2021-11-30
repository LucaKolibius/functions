clear
folderpath = '\\analyse4.psy.gla.ac.uk\project0309\Luca\official code'; % wherever you have downloaded this folder
load([folderpath, '\data\allSpks.mat'], 'allSpks');

figure(1); clf;
for tarNum = [10 15] % number of TC
    allIU             = [allSpks.iu];
    allIU(allIU == 2 | allIU == 1) = 1;
    nperm             = 10000;
    % tarNum            = 9; % [3 IU / 9 TC] [3 iu / 15 TC]
    
    % How many IU are expected to be drawn from tarNum = 16 drawings under the H0?
    for perm = 1:nperm
        
        allIU         = allIU(randperm(size(allIU,2)));
        permNum(perm) = sum(allIU(1:tarNum));
        
    end
    
    numIU = 1;
    pVal = mean(permNum >= numIU)
        
    switch tarNum
%         case 9
%             subplot(211)
%             histogram(permNum, 0:1:10, 'normalization', 'probability')
% %             endP =  sum(permNum>=3) / 10000;
%             endP = mean(permNum>4);
%             title(sprintf('Number of IU in a randomly drawn sample of %d | p = %.2f for three IU or more', tarNum,endP))
%             ylabel('Probability')
%             xlabel('Number of IU')
%             ylim([0 0.4])
        case 16
            subplot(212)
            histogram(permNum, 0:1:10, 'normalization', 'probability')
%             endP =  sum(permNum>=3) / 10000;
endP = mean(permNum>4);
            title(sprintf('Number of IU in a randomly drawn sample of %d | p = %.2f for three IU or more', tarNum,endP))
            ylabel('Probability')
            xlabel('Number of IU')
            ylim([0 0.4])
    end
end


