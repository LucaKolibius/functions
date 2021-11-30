% tic
% V = a - b;
% D = sqrt(V * V');
% toc

T      = [];
T(1,:) = randn(1,1990)*0.2;
T(2,:) = randn(1,1990)*0.2;

% a = [(T(1,:)*5)-1; (T(2,:)*5)];
% b = T;
% b = [];
% c = [T(1,:)+1; T(2,:)];

distShift = [-0.5 +0.25];
a = [(T(1,:)*5)+distShift(1); (T(2,:)*5)];
b = T;
b = [];
c = [T(1,:)+distShift(2); T(2,:)];


allP = [a b c];

tic
count = 1;
D = [];
for c1 = 1:length(allP)
    for c2 = 1:length(allP)
        
        if c1==c2
            continue
        end
        
        V = allP(:,c1) - allP(:,c2);
        D(count) = sqrt(V'* V);
        count = count + 1;
    end
end
toc

figure(1); clf;
subplot(211); hold on;
scatter(a(1,:), a(2,:))
% scatter(b(1,:), b(2,:))
scatter(c(1,:), c(2,:))

subplot(212)
histogram(D)

tic
[~, clPos] = kmeans(allP', 2);
toc

subplot(211); hold on;
scatter(distShift(1), 0, 100, 'k', 'filled');
scatter(distShift(2), 0, 100, 'k', 'filled');

% plot(clPos(1,1), clPos(1,2),'-ro', 'MarkerSize',20)
scatter(clPos(1,1), clPos(1,2), 60, 'r', 'filled');
scatter(clPos(2,1), clPos(2,2), 60, 'r', 'filled');


findpeaks % minpeakwith? at least minpeakprominence!