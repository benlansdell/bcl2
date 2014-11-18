%simple_bistable.m
%Script to demonstrate the bistability phenomenona caused by auto-activation
%Without bistability:

%Set parameters
t = 0:0.1:10;
ic = 1;
K = 4;
df = @(t, x) K-x*(1+K);

%Numerical solution:
[s, f] = ode15s(df, t, ic);

%Exact solution:
exact = @(t) (K + exp(-(1+K)*t))/(1+K);

%Plot exact and numerical solution to check they are the same...
plot(t, f, t, exact(t));

%Save image
filename = './images/demos/simple_bistable-1.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%With bistability
%Set parameters
t = 0:0.01:9.9;
ic = 1; k1 = 1;
k2 = 2; k3 = 3;
B0 = 10;
sol = [];

for T0 = 100
	%Numerical solution
	K2 = k2/k1/T0;
	K3 = B0*k3/k1/T0;
	df = @(t, x) K2 - x*(1 + K3 + K2) + K3*x.*x;
	[s, f] = ode15s(df, t, ic);
	sol = [sol; f'];
end

%Exact solution:
%c0 = 1/K3 + 1 + K2/K3
%d0 = sqrt(c0*c0/4 - K2/K3)
%A = (1- d0 - c0/2)/(1 - c0/2 + d0)
%exact = @(t) d0*(exp(-2*d0*K3*t)+A)./(exp(-2*d0*K3*t)-A) + c0/2;

%Plot exact and numerical solution to check they are the same...
figure
plot(t, sol);
xlabel('time (s)');
ylabel('\beta');

%Save image
filename = './images/demos/simple_bistable-2.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

%Compute bifurcation diagram for bistable case...
r = [];
for T0 = 0:20
	K2 = k2/k1/(T0+.1);
	K3 = B0*k3/k1/(T0+.1);
	%Solve polynomial
	r = [r roots([K3 (-1 - K3 - K2) K2])];
end
r
plot(0:20, r)
xlabel('T0');
ylabel('\beta');
filename = './images/demos/simple_bistable-3.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');

