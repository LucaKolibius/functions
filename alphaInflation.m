clear
thresh = [];
nperm = 10000;
numSU = 300;
tic

for wdh = 1:1000
    toc
    clearvars -except nperm numSU plvl wdh empDat
    trls = 50;
    ewpPerm = zeros(50,nperm,numSU);
    thresh  = [];
    disp(wdh);
    for su  = 1:numSU
        %disp(su)
        enc = randn(trls,1);
        ret = randn(trls,1);
        ewp = enc .* ret;
        
        for perm = 1 : nperm
            encPerm            = enc(randperm(trls));
            retPerm            = ret(randperm(trls));
            ewpPerm(:,perm,su) = encPerm .* retPerm;
        end
        
        tmp       = ewpPerm(:,:,su);
        thresh(su) = prctile(tmp(:), 99);
        
        if sum(ewp >= thresh(su)) >= 2
            empSig(su) = 1;
        else
            empSig(su) = 0;
        end
       % toc
    end % all numSU

    % second level permutation
    for perm = 1:nperm
        temp = [];
        for su = 1:numSU
            
            % take a random permutation
            permNum  = randperm(nperm,1);
            thisPerm = ewpPerm(:,permNum,su);
                
            if sum(thisPerm >= thresh(su)) >= 2
                temp(su) = 1;
            else
                temp(su) = 0;
            end
               
        end
        
        numSig(perm) = sum(temp); % number of "index units" under the null
        
    end
    empDat (wdh) = sum(empSig);
%     permDat(:, wdt) = prctile(numSig,95,2);
    
    plvl(wdh) = 1- sum(empDat(wdh) >= numSig) / nperm;
end % on to the next repetition
save('alphaOP3.mat', 'plvl', 'empDat');
