%Script to plot steady states as computed by reversible_steady_states.mpl in Maple...
%Load matrix files
load scripts/reversible_steady_states-tBid-1.dat;
load scripts/reversible_steady_states-tBIM-1.dat;

tBid = reversible_steady_states_tBid_1;
tBIM = reversible_steady_states_tBIM_1;
len = length(tBIM);

%Plot b
plot(1:len, tBid(1,:), 1:len, tBIM(1,:))
xlabel('T_0');
ylabel('\beta');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-1-1.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot t
plot(1:len, tBid(2,:), 1:len, tBIM(2,:))
xlabel('T_0');
ylabel('\tau');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-1-2.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot m
plot(1:len, tBid(3,:), 1:len, tBIM(3,:))
xlabel('T_0');
ylabel('\mu');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-1-3.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot d
plot(1:len, tBid(4,:), 1:len, tBIM(4,:))
xlabel('T_0');
ylabel('\delta');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-1-4.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot b_hat
plot(1:len, tBid(5,:), 1:len, tBIM(5,:))
xlabel('T_0');
ylabel('\beta bar');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-1-5.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot t_hat
plot(1:len, tBid(6,:), 1:len, tBIM(6,:))
xlabel('T_0');
ylabel('\tau bar');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-1-6.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot m_hat
plot(1:len, tBid(7,:), 1:len, tBIM(7,:))
xlabel('T_0');
ylabel('\mu bar');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-1-7.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

