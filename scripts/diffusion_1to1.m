%Plot solutions to simple 1:1 binding model which includes diffusive effects.
%Use dimensionless variables

%Dimesional parameters
%Bak:Mcl-1 Biacore on and off rates
k_on = 7.98e-4;
k_off = 1.37e-3;
%Bak
A_0 = 10;
%Mcl-1
B_0 = 20;

%Dimensionless parameters
sigma = A_0/B_0
K = k_off/k_on/B_0

equil = C_0/2 - D_0
t = 0:0.1:10;
ic = [0];
sol = [];

%Try different D values:
for D = [0 10^(-7.5) 10^(-8) 10^(-8.5) 10^(-9)]
	d = A_0*k_on*5/6e6/D;
	v1 = @(x, sigma, K, d) (1-x(1)).*(1-sigma*x(1))./(1+d*(1-x(1))) - K*x(1)./(1+d*(1-x(1)));
	df = @(t, x) [v1(x, sigma, K, d)];
	[s, y] = ode15s(df, t, ic);
	sol = [sol; y'];
end

%Compute the others from conservations
sol_a = 1 - sol;
sol_b = 1 - sigma*sol;

%Compute the exact solution also...
C_0 = 1/sigma + K/sigma + 1;
D_0 = sqrt(C_0^2/4 - 1/sigma);
A = (C_0 + 2*D_0)/(C_0 - 2*D_0);
exact = @(t, sigma, C_0, D_0) D_0*(exp(-2*D_0*t*sigma) + A)./(exp(-2*D_0*t*sigma) - A) + C_0/2;

figure
prev_ColorOrder = get(0, 'DefaultAxesColorOrder');
my_ColorOrder = [0.00000   0.00000   0.00000;
			     0.00000   0.00000   0.25000;
			     0.00000   0.00000   0.50000;
  			     0.00000   0.00000   0.75000;
  				 0.00000   0.00000   1.00000;];
set(0,'DefaultAxesColorOrder',my_ColorOrder);

plot(t, sol)
%line(t, exact(t, sigma, C_0, D_0), 'Color', 'black');
line(t, equil, 'Color', 'black', 'LineStyle', '--');
xlabel('\tau');
ylabel('concentration');
ylim([0 1]);
legend('d = 0 (reaction-limited)', 'd = 0.2103', 'd = 0.6650', 'd = 2.1029', 'd = 6.6500', 'Location', 'SouthEast');

filename = './images/demos/diffusion_1to1.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');
set(0, 'DefaultAxesColorOrder', prev_ColorOrder);