%Script to compare the full michaelis-menton model of a catalysed reaction to the simple one:
%Full:
% A + B <-> AB --> C + B

a_conc = [];
c_conc = [];

%Set parameters
t = 0:0.1:100;
ic = [5; 2; 0; 0];

for i = 1:6
	k = [1e-2 1e-1; 10^(i/2-2) 0];
	v1 = @(x, k, ic, t) k(1,1)*x(1)*x(2) - k(1,2)*x(3);
	v2 = @(x, k, ic, t) k(2,1)*x(3);
	df = @(t, x) [-v1(x, k, ic, t); -v1(x, k, ic, t) + v2(x, k, ic, t); v1(x, k, ic, t) - v2(x, k, ic, t); v2(x, k, ic, t)];

	%Numerical solution
	[s, f] = ode15s(df, t, ic);
	a_conc = [a_conc; f(:,1)'];
	c_conc = [c_conc; f(:,4)'];
end

len = length(t);
prev_ColorOrder = get(0, 'DefaultAxesColorOrder');
my_ColorOrder = [0.00000   0.00000   0.00000;
			     0.00000   0.00000   0.20000;
			     0.00000   0.00000   0.40000;
  			     0.00000   0.00000   0.60000;
  				 0.00000   0.00000   0.80000;
  				 0.00000   0.00000   1.00000;
  				 0.00000   0.00000   0.00000;
			     0.20000   0.00000   0.00000;
			     0.40000   0.00000   0.00000;
  			     0.60000   0.00000   0.00000;
  			     0.80000   0.00000   0.00000;
  				 1.00000   0.00000   0.00000;];
set(0,'DefaultAxesColorOrder',my_ColorOrder);

figure
h=plot(t/60, a_conc, t/60, c_conc);
order = [6 12];
xlabel('time (min)');
ylabel('concentration (nM)');
ylim([0 6]);
%legend(h(order), {'A', 'C'});
text(.1, 5, 'A');
text(.1, 1, 'C');

filename = './images/demos/simple_michaelis_menton.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');
set(0, 'DefaultAxesColorOrder', prev_ColorOrder);