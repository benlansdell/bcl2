%Plot solution to simple 1:1 binding model as described in Example 3.1

t = 0:0.1:100;
k = [1e-2 1e-1];
ic = [5; 2; 0];
v1 = @(x, k, ic, t) k(1,1)*x(1)*x(2) - k(1,2)*x(3);
df = @(t, x) [-v1(x, k, ic, t); -v1(x, k, ic, t); v1(x, k, ic, t)];
%Numerical solution
[s, f] = ode15s(df, t, ic);

%`Exact' solution 1
f1 = @(k, ic, t) (ic(1) - ic(2))*ic(2)*exp(-k(1,1)*(ic(1)-ic(2))*t)./(ic(1) - ic(2)*exp(-k(1,1)*(ic(1)-ic(2))*t));

%`Exact' solution 2
c0 = @(k, ic) ic(1) + ic(2) + k(1,2)/k(1,1);
d0 = @(k, ic) sqrt(c0(k, ic)*c0(k, ic)/4 - ic(1)*ic(2));
A = @(k, ic) (d0(k, ic) + c0(k, ic)/2)/(c0(k, ic)/2-d0(k, ic));
f2 = @(k, ic, t) d0(k, ic)*(exp(-2*d0(k, ic)*k(1,1)*t)+A(k, ic))./(exp(-2*d0(k, ic)*k(1,1)*t)-A(k, ic)) + c0(k, ic)/2;

%`Exact' solution 3
f3 = @(k, ic, t) (d0(k, ic)+.5*c0(k, ic))*(1+(2*d0(k, ic))./(exp(-2*d0(k, ic)*k(1,1)*t)+1-c0(k, ic)));

len = length(t);
figure
plot(t, f)
%line(t, f1(k, ic, t));
line(t, f2(k, ic, t));
%line(t, f3(k, ic, t));
xlabel('time (s)');
ylabel('concentration (nM)');
ylim([0 5.5]);
%legend('A', 'B', 'C', 'A1', 'C2', 'C3_longterm');
legend('A', 'B', 'C');

filename = './images/demos/example3_1.eps';
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [6 4]);
set(gcf, 'paperposition', [0 0 6 4]);
print(gcf,  filename, '-depsc');