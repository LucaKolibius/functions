function spikeNum = insertSpikeNumber(trig, spks, tw)
spikeNum = zeros(length(trig), 1);

for trl = 1 : length(trig)
    
    spkTrl          = spks >= trig(trl)+tw(1) & spks <= trig(trl)+tw(2);
    spikeNum(trl,1) = sum(spkTrl);
    
end % end of function