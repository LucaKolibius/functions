doAnal = true(625,1);

allBids = {allSpks.bidsID};
allSesh = {allSpks.sesh};
output  = struct;
count = 1;
for su = 1:625
    
    if doAnal(su) == false
        continue
    end
    bidsID                = allSpks(su).bidsID;
    sesh                  = allSpks(su).sesh;
    
    sameSesh              = strcmp(allBids, bidsID) & strcmp(allSesh, sesh);
    doAnal(sameSesh)      = false;
    
    reinstSesh            = vertcat(allSpks(sameSesh).reinstTrl);
    reinstSesh            = reinstSesh(:,allSpks(su).hitsIdx);
    
    output(count).numEv   = size(reinstSesh,2);
    output(count).numSu   = size(reinstSesh,1);
    output(count).numESN  = sum(reinstSesh,1);
    output(count).propESN = mean(reinstSesh,1);
    output(count).reinEv  = any(reinstSesh,1);
    count = count + 1;
  
end

numEv   = [output(:).numEv];
numSu   = [output(:).numSu];
numESN  = [output(:).numESN];
propESN = [output(:).propESN];
reinstEv = sum([output(:).reinEv]);

numESN(numESN==0)   = [];
propESN(propESN==0) = [];

