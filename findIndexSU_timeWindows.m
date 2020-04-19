function [timeWindowENC, timeWindowRET, lockingENC, lockingRET, onlyCues] = findIndexSU_timeWindows(encWindow, retWindow)
onlyCue(1:27) = 0;

timeWindows{1}   = [-1  0];
lockings{1}      = [ 1  1];
onlyCue(1)       = 1;

timeWindows{2}   = [-1  1];
lockings{2}      = [ 1  1];
onlyCue(2)       = 1;

timeWindows{3}   = [-1  2];
lockings{3}      = [ 1  1];
onlyCue(3)       = 1;

timeWindows{4}   = [-1  -2];
lockings{4}      = [ 1   3];

timeWindows{5}   = [-1 -1];
lockings{5}      = [ 1  3];

timeWindows{6}   = [-1  0];
lockings{6}      = [ 1  3];

timeWindows{7}   = [-1  1];
lockings{7}      = [ 1  3];

timeWindows{8}   = [ 0  1];
lockings{8}      = [ 1  1];
onlyCue(8)       = 1;

timeWindows{9}   = [ 0  2];
lockings{9}      = [ 1  1];
onlyCue(9)       = 1;

timeWindows{10}  = [ 0 -2];
lockings{10}     = [ 1  3];

timeWindows{11}  = [ 0 -1];
lockings{11}     = [ 1  3];

timeWindows{12}  = [ 0  0];
lockings{12}     = [ 1  3];

timeWindows{13}  = [ 0  1];
lockings{13}     = [ 1  3];

timeWindows{14}  = [ 1  2 ];
lockings{14}     = [ 1  1 ];
onlyCue(14)      = 1;

timeWindows{15}  = [ 1 -2];
lockings{15}     = [ 1  3];

timeWindows{16}  = [ 1 -1];
lockings{16}     = [ 1  3];

timeWindows{17}  = [ 1  0];
lockings{17}     = [ 1  3];

timeWindows{18}  = [ 1  1];
lockings{18}     = [ 1  3];

timeWindows{19}  = [ 2 -2];
lockings{19}     = [ 1  3];

timeWindows{20}  = [ 2 -1];
lockings{20}     = [ 1  3];

timeWindows{21}  = [ 2  0];
lockings{21}     = [ 1  3];

timeWindows{22}  = [ 2  1];
lockings{22}     = [ 1  3];

% timeWindows{23}  = [-2 -1];
% lockings{23}  = [ 3  3];
% 
% timeWindows{24}  = [-2  0];
% lockings{24}  = [ 3  3];
% 
% timeWindows{25}  = [-2  1];
% lockings{25}  = [ 3  3];
% 
% timeWindows{26}  = [-1  0];
% lockings{26}  = [ 3  3];
% 
% timeWindows{27}  = [-1  1];
% lockings{27}  = [ 3  3];
%%
timeWindowENC =  timeWindows{encWindow};
lockingENC    =  lockings{encWindow};
timeWindowRET =  timeWindows{retWindow};
lockingRET    =  lockings{retWindow};
onlyCues      =  onlyCue(encWindow) + onlyCue(retWindow);

% only if both timewindows (for encoding and retrieval) only include the
% cue, then onlyCue is considered (no other stimuli but the cue were shown,
% making single trial index neurons, as opposed to double trial index
% neurons, possible)
if onlyCues == 2
    onlyCues = 1;
elseif onlyCues < 2
    onlyCues = 0;
end

end % end of function        
        