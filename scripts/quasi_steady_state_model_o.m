%Find steady state concentrations when assume that there is no inactive Bak, nor does any Bak
%dimerisation occur. Hence, we have Bak*:Mcl-1 complex and tBid:Mcl-1 complex non-trival concentrations.

%From derivations we solve the following:

%Starting concentrations
B0 = 10;
M0 = 20;
T0 = 300;

%Non-dimensional
alpha = T0/M0;
gamma = M0/B0;

%Biacore numbers (for tBid)
k2a = 1.37e-3;		%Bak:Mcl-1 assoc (1/nM/s)
k2d = 7.98e-4;		%Bak:Mcl-1 dissoc (1/s)
k3a = 8.10e-3;	    %tBid:Mcl-1 assoc (1/nM/s)
k3d = 7.52e-2;	    %tBid:Mcl-1 dissoc (1/s)

%Non-dimensional
K2D = k2d/k2a/B0;
K3D = k3d/k3a/M0;

%Our calculations show
a_1 = -alpha*K2D + alpha*gamma*K3D;
a_2 = 2*K2D*alpha + K2D + K3D*K2D + K3D - K3D*gamma - K3D*alpha*gamma - gamma*K3D^2;
a_3 = -2*K2D - alpha*K2D - K2D*K3D - K3D + gamma*K3D;
a_4 = K2D;

%Solve it
p = [a_1 a_2 a_3 a_4];
r = roots(p);
r = r(r<=1&r>=0);

if (length(r) ~= 1)
	throw(MException('RootErr:TooManyTooFew', 'Too many or too few valid steady state concentrations'));
end

%And rescale
tau_hat = r(1);
mu = K3D*tau_hat/(1-tau_hat);
mu_hat = 1-mu-alpha*tau_hat;
tau = 1-tau_hat;
beta_hat = 1-gamma*mu_hat;

T = T0*tau
T_hat = T0*tau_hat
B_hat = B0*beta_hat
M = M0*mu
M_hat = M0*mu_hat