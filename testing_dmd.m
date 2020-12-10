%%
clear
load('X:\Luca\data\microLFP\sub-0007_S1_onlyMicroLFP_RAW_1000DS_SPKINT.mat', 'data')

LFP     = data.trial{1}(:,1:10000);
dt      = data.time{1}(2) - data.time{1}(1);
origRow = size(LFP,1);

diffLFP = zeros(438,50);
counter = 1;

lfpIdxAll = round(unifrnd(repmat(2000,1,20),8000));

bestPred = 99999;
for winLen = 40:20:1000
    disp(winLen)
    for nstacks = 10:5:50
        tempDiff = [];
        for iter = 1:20
            lfpIdx = lfpIdxAll(iter);
            
            try
                LFPwin             = LFP(:,lfpIdx-winLen:lfpIdx-1);
                [~, ~, ~, ~, Xdmd] = DMDfull(LFPwin, {'dt', 1, 'nstacks', nstacks});
                
                %% unstack
                unstacked = Xdmd(1:origRow,:);
                for stack = 2:nstacks
                    unstacked = [unstacked, Xdmd(origRow*(stack-1)+1:origRow*stack,:)];
                end
                
                reconLFP           = real(unstacked);
                reconLFPpred       = reconLFP(:,winLen+1:winLen+50);
                LFPreal            = LFP(:,lfpIdx:lfpIdx+49);
                
                
                tempDiff(iter,:)   = mean(abs(LFPreal-reconLFPpred),1);
                
                
                if sum(tempDiff(iter,:),2) < bestPred
                    saveLFP.recon = reconLFPpred;
                    saveLFP.real  = LFPreal;
                    bestPred = sum(tempDiff(iter,:),2);
                end
            catch
                continue
            end
            
               
                
            
        end
        if isempty(tempDiff)
            continue
        end
        diffLFP(counter,:) = mean(tempDiff,1);
        output(1,counter)  = winLen;
        output(2,counter)  = nstacks;
        counter            = counter+1;

    end
end
    
    
    
    figure(1); clf;
    hold on
    plot(reconLFPpred(5,:), 'k');
    plot(LFPreal(5,:), 'r');
    
    figure(1); clf;
    hold on;
    plot(reconLFP(1,1:winLen), 'r');
    plot(LFPwin(1,:), 'k');
    
    
    % plot(imagP(2,:), 'b');
    % subplot(311);
    % normLFP = LFP / max(LFP(:));
    % normLFP = normLFP + repmat([1:32]', 1,1000);
    % plot(normLFP', 'color', 'k');
    %
    % subplot(312);
    % normLFP = realP / max(realP(:));
    % normLFP = normLFP + repmat([1:32]', 1,999);
    % plot(normLFP', 'color', 'k');
    
    
    % , {'r', [1e32]}, {'nsttacks', 1}
    
    %   {'dt', [1]}
    %   {'r', [1e32]}        truncate to rank-r
    %   {'nstacks', 1}       number of stacks of the raw data
    %% stacking
    nstacks = 30; % number of stacks
    
    % construct the augmented, shift-stacked data matrices
    Xstack = [];
    for st = 1 : nstacks
        Xstack = [Xstack; LFP(:, st:end-nstacks+st)];
    end
    
    X1 = Xstack(:, 1:end-1);
    X2 = Xstack(:, 2:end);
    m = size(X1,1);
    
    % SVD and truncate to first r modes
    r = 100;
    [U, S, V] = svd(X1, 'econ');
    U_r = U(:, 1:r);
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);
    
    Atilde = U_r' * X2 * V_r / S_r;
    [W_r, D] = eig(Atilde);
    Phi = X2 * V_r / S_r * W_r;
    
    lambda = diag(D);
    omega = log(lambda)/dt;
    
    x1 = X1(:, 1);
    b = Phi\x1;
    
    %%
    
    t = (1:m)/dt;
    time_dynamics = [];
    for iter = 1:m
        time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
    end
    Xdmd = Phi * time_dynamics;
    Xdmd = real(Xdmd);
    
    test = Phi*lambda*exp(961)*b
    
    
    % b = Phi / Atilde; ???
    predict = zeros(1000,1000);
    for tmp = 1:1000
        for dmd = 1 : r
            predict(:,tmp) = predict(:,tmp)+(Phi(dmd,:)*exp(W_r(:,dmd)*tmp)*b(dmd,:))';
        end
    end
    
    
    Phi*lambda(r)*b
    
    for ii = 1:3
        xdmd(ii,:) = Phi(:,ii)*omega(ii)*b(ii);
    end
    
    xdmd = sum(xdmd,1);
    xdmd = abs(xdmd);
    xdmd = xdmd - mean(xdmd,2);
    
    pw = abs(fft(xdmd));
    hz = linspace(0,20,1000);
    
    figure(1);clf;
    subplot(211);
    plot(xdmd(1:100))
    
    subplot(212);
    plot(hz,pw);
    
