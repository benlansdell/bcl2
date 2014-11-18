#Script to reduce the dimensionless reversible Bcl2 model 
#to only beta, tau, mu, delta, gamma variables
#Allowing eigenvalues for stability matrix to be calculated...

with(Statistics);
with(ArrayTools);
with(LinearAlgebra);
with(VectorCalculus);
with(StringTools);

if (Compare(activator, 'tBIM')) then	
	#Values from fitted parameters (tBIM case)
	k1a := 10^(-5.83);
	#k1d := 10^(-6);
	k2a := 10^(-3.46);
	k2d := 10^(-3.21);
	k3a := 10^(-5.24);
	k3d := 10^(-6);
	k4a := 10^(-3.57);
	k4d := 10^(-4);
	k5a := 10^(-4.11);
elif (Compare(activator, 'tBid')) then
	#Values from fitted parameters (tBid case)
	k1a := 10^(-3.51);
	#k1d := 10^(-6);
	k2a := 10^(-3.46);
	k2d := 10^(-3.21);
	k3a := 10^(-4.84);
	k3d := 10^(-3.6);
	k4a := 10^(-3.57);
	k4d := 10^(-4);
	k5a := 10^(-4.11);
else
	return(FAIL);
end if:	

M0 := 20;
B0 := 10;
T0 := 30;
g := M0/B0;
a := T0/M0;

K1D := k1d/k1a/T0;
K2A := k2a*B0/k1a/T0;
K2D := k2d/k1a/T0;
K3A := k3a*M0/k1a/T0;
K3D := k3d/k1a/T0;
K4A := k4a*B0/k1a/T0;
K4D := k4d/k1a/T0;
K5A := k5a*B0/k1a/T0;

b_hat := K2D*m_hat/K2A/m;
m     := K3D*(1-t)/K3A/t;
d     := K4A*b_hat*b_hat/K4D;
b     := K1D*b_hat/(t+K5A*b_hat);
t_hat := 1-t;
m_hat := 1-m-a+a*t;
f := b_hat + b + d + g*m_hat;

biol := (x,y) -> if (Im(eval(x, t=y)) <> 0 or eval(x,t=y)<0 or eval(x,t=y)>1) then false else true fi;

tss := []:
tss_ := [solve(f=1, t)]:
for ss in tss_ do
	#For each ss check it is biologically meaningful
	#That is, every concentration for every species must be real and between 0 and 1
	if (biol(t, ss) and biol(b, ss) and biol(m, ss) and biol(d, ss) and biol(b_hat, ss) and biol(t_hat, ss) and biol(m_hat, ss)) then
		tss := [op(tss), ss]:
	end if:
end do:

b_ss := Vector[column](map(x->eval(b,t=x), tss)));
t_ss := Vector[column](tss)):
m_ss := Vector[column](map(x->eval(m,t=x), tss))):
d_ss := Vector[column](map(x->eval(d,t=x), tss))):
b_hat_ss := Vector[column](map(x->eval(b_hat,t=x), tss))):
t_hat_ss := Vector[column](map(x->eval(t_hat,t=x), tss))):
m_hat_ss := Vector[column](map(x->eval(m_hat,t=x), tss))):

S := Matrix([[-1,0,0,0,-1,0],[0,0,-1,0,0,0],[0,-1,-a,0,0,0],[0,0,0,1,0,0],[0,0,0,0,0,-1]]);
J := Vector([b*t - K1D*b_hat,K2A*m*b_hat-K2D*m_hat,K3A*t*m-K3D*t_hat,K4A*b_hat*b_hat-K4D*d,K5A*b*b_hat,K6A*d*g]);
j := Jacobian(J, [b,t,m,d,g]=[bs,ts,ms,ds,gs]);
A := Multiply(S,j);
evals := Eigenvalues(A);
