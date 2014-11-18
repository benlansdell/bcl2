%Script to plot steady states as computed by reversible_steady_states.mpl in Maple...
%Load matrix files
load scripts/reversible_steady_states-tBid-2.dat;
load scripts/reversible_steady_states-tBIM-2.dat;

tBid = reversible_steady_states_tBid_2;
tBIM = reversible_steady_states_tBIM_2;
len = length(tBIM);

%Plot b
plot((1:len)/4-7, tBid(1,:), (1:len)/4-7, tBIM(1,:))
xlabel('k4d');
ylabel('\beta');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-2-1.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot t
plot((1:len)/4-7, tBid(2,:), (1:len)/4-7, tBIM(2,:))
xlabel('k4d');
ylabel('\tau');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-2-2.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot m
plot((1:len)/4-7, tBid(3,:), (1:len)/4-7, tBIM(3,:))
xlabel('k4d');
ylabel('\mu');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-2-3.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot d
plot((1:len)/4-7, tBid(4,:), (1:len)/4-7, tBIM(4,:))
xlabel('k4d');
ylabel('\delta');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-2-4.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot b_hat
plot((1:len)/4-7, tBid(5,:), (1:len)/4-7, tBIM(5,:))
xlabel('k4d');
ylabel('\beta bar');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-2-5.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot t_hat
plot((1:len)/4-7, tBid(6,:), (1:len)/4-7, tBIM(6,:))
xlabel('k4d');
ylabel('\tau bar');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-2-6.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Plot m_hat
plot((1:len)/4-7, tBid(7,:), (1:len)/4-7, tBIM(7,:))
xlabel('k4d');
ylabel('\mu bar');
legend('tBid', 'tBIM');
filename = './images/demos/reversible_steady_states-2-7.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

