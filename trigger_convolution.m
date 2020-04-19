function [ xy_d, out ] = trigger_convolution( trg_1, trg_2 )
% Convolve y across x, where x is longer than y. If out is 1, then trg_1 is
% longer, if out is -1, then trg_2 is longer.
if(length(trg_1)>=length(trg_2)); x = trg_1; y = trg_2; out = 1; else; x = trg_2; y = trg_1; out = -1; end

Ly = length(y); Lx = length(x);
Ld = length(Ly/2 : Lx - Ly/2);
xy_d = zeros(1, Ld);
y = y - min(y);

for i = Ly/2 : Lx - Ly/2
    x_t = x(i - Ly/2 +1 : i + Ly/2); 
    x_t = x_t - min(x_t);
    xy_d(i - Ly/2 + 1) = mean(abs(x_t - y)); 
end

if min(xy_d) < 0.05
    xy_d = find(xy_d == min(xy_d));
else
    error('TTL-Logfile Match not good enough. Investigate.')
end
end

