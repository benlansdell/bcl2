%Script to plot steady states as computed by reversible_steady_states.mpl in Maple...
%Load matrix files
load scripts/reversible_steady_states-tBid.dat;
load scripts/reversible_steady_states-tBIM.dat;

tBid = reversible_steady_states_tBid;
tBIM = reversible_steady_states_tBIM;
len = length(tBIM);

%Plot b
plot(1:len, tBid(1,:), 1:len, tBIM(1,:))
xlabel('T0');
ylabel('\beta');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-1.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot t
plot(1:len, tBid(2,:), 1:len, tBIM(2,:))
xlabel('T0');
ylabel('\tau');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-2.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot m
plot(1:len, tBid(3,:), 1:len, tBIM(3,:))
xlabel('T0');
ylabel('\mu');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-3.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot d
plot(1:len, tBid(4,:), 1:len, tBIM(4,:))
xlabel('T0');
ylabel('\delta');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-4.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot b_hat
plot(1:len, tBid(5,:), 1:len, tBIM(5,:))
xlabel('T0');
ylabel('\beta bar');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-5.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot t_hat
plot(1:len, tBid(6,:), 1:len, tBIM(6,:))
xlabel('T0');
ylabel('\tau bar');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-6.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot m_hat
plot(1:len, tBid(7,:), 1:len, tBIM(7,:))
xlabel('T0');
ylabel('\mu bar');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-7.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

