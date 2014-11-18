[k, ic] = load_model('tBim', 'reversible_exp4_BIAcore');

B0 = ic(1);
T0 = ic(2);
M0 = ic(3);
C0 = ic(5);

%Provided that parameters are loaded from load_model() and adjusted by fit_exp4_*() then
%we know the parameters are non-zero and that cytochrome c release is therefore
%'inevitable'. Thus...
c = 0;
ch = C0;

%Non-dimensionalise the system
alpha = T0/M0;
lambda = M0/B0;
K1D = k(1,2)/k(1,1)/T0;
K2A = k(2,1)*B0/k(1,1)/T0;
K2D = k(2,2)/k(1,1)/T0;
K3A = k(3,1)*M0/k(1,1)/T0;
K3D = k(3,2)/k(1,1)/T0;
K4A = 2*k(4,1)*B0/k(1,1)/T0;
K4D = k(4,2)/k(1,1)/T0;
K5A = k(5,1)*B0/k(1,1)/T0;
K6A = k(6,1)*B0/k(1,1)/T0;
K6D = k(6,2)/k(1,1)/T0;
K7A = k(7,1)*B0/k(1,1)/T0/4;

%Solve dimensionless system first...very slow
syms b t m d bh th mh b4;
S = solve(t+th-1, m+mh+alpha*th-1, b+bh+d+b4+lambda*mh-1, K2A*bh*m-K2D*mh, K3A*t*m-K3D*th, K6A*d*d-K6D*b4, K4A*bh*bh-K4D*d, b*t-K1D*bh+K5A*b*bh);
solution = physical_sol(S);

%Try simplifying the solution first and then creating a polynomial to solve numerically...maybe quicker
digits(5);
a1 = subs(b*t + K1D*bh - K5A*b*bh, {b}, {1-bh-d-b4+lambda*mh})
a2 = subs(a1, {b4}, {K6A*d*d/K6D})
a3 = subs(a2, {d}, {K4A*bh*bh/K4D})
a4 = subs(a3, {mh}, {K2A*bh*m/K2D})
a5 = subs(a4, {m}, {K3D*(1-t)/K3A/t})
a6 = subs(a5, {bh}, {K2D*K3A*t*(1 - K3D*(1-t)/t/K3A - alpha + alpha*t)/K2A/K3D/(1-t)})
a7 = vpa(expand(a6), 4);
