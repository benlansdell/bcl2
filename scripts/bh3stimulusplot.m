%Plot solution to simple 1:1 binding model as described in Example 3.1
t = 0:0.1:3600;
window = 60;
%New tBim added to system...

%Exp 4
%T_added_prime = @(t) (0.7*(t>600-window).*(t<=600) + 2*(t>1200-window).*(t<=1200) + 7*(t>1800-window).*(t<=1800) + 20*(t>2400-window).*(t<=2400) + 70*(t>3000-window).*(t<=3000))/window;
%T_added = @(t) 0.3 + 0.7*min(max(0, (t-600+window)/window), 1) + 2*min(max(0, (t-1200+window)/window), 1) + 7*min(max(0, (t-1800+window)/window), 1) + 20*min(max(0, (t-2400+window)/window), 1) + 70*min(max(0, (t-3000+window)/window), 1);

%Exp 6
T_added_prime = @(t) (7*(t>600-window).*(t<=600) + 20*(t>1200-window).*(t<=1200) + 70*(t>1800-window).*(t<=1800) + 200*(t>2400-window).*(t<=2400))/window;
T_added = @(t) 3 + 7*min(max(0, (t-600+window)/window), 1) + 20*min(max(0, (t-1200+window)/window), 1) + 70*min(max(0, (t-1800+window)/window), 1) + 200*min(max(0, (t-2400+window)/window), 1);

len = length(t);
figure
plot(t/60, T_added(t))
xlabel('time (min)');
ylabel('concentration (nM)');
ylim([0 310]);
saveplot(gcf, './images/demos/bh3stimulusplot1.eps');
figure
plot(t/60, T_added_prime(t))
xlabel('time (min)');
ylabel('concentration/second (nM s^{-1})');
ylim([0 3.5]);
saveplot(gcf, './images/demos/bh3stimulusplot2.eps');