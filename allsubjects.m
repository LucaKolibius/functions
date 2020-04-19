% %%
% clear all; subjID = 'P01ERL_S1';
% disp(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
% RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P01ERL_S2';
% disp(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
% RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

% %%
% clear all
% subjID = 'P02_S1';
% disp(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
% RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P02_S2';
% disp(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
% RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P02_S3';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
% RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P02_S4';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
% RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

% clear all; subjID = 'P02_S5';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
% % RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P02_S6';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
% RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P02_S7';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
% RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

%% P02ERL IS FUCKING WEIRD! TTLS AND TRIGGERS DONT LINE UP / GEORGES CROSSCORRELATION?
clear all; subjID = 'P02ERL_S1';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P02ERL_S2';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P03ERL_S1';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P03ERL_S2';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P03ERL_S3';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P03ERL_S4';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

%% ONLY 1 CL
% MB: left anterior hipp; MC: left post hip // no idea which MACRO
clear all; subjID = 'P04ERL_S1';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

% 3 CL
clear all; subjID = 'P04ERL_S2';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P04ERL_S3';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

%% no localization available
% S1 all CL fine
clear all; subjID = 'P04_S1';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P04_S2';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P04_S3';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

%% 
clear all; subjID = 'P05ERL_S1';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P05ERL_S2';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P05ERL_S3';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P05ERL_S4';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

%%
clear all; subjID = 'P05_S1';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P05_S2';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P05_S3';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

clear all; subjID = 'P05_S4';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
    
clear all; subjID = 'P05_S5';
disp(subjID);
rename_clNames(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID);  
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 

% %%
% clear all; subjID = 'P07_S1';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P07_S2';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P07_S3';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% %%
% clear all; subjID = 'P08_S1';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P08_S2';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P08_S3';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P08_S4';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% %%
% clear all; subjID = 'P09_S1';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P09_S2';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 
% 
% clear all; subjID = 'P09_S3';
% disp(subjID);
% rename_clNames(subjID);
% [encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID); 