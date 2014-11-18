#Script to compute the steady state(s?) of the reversible Bcl2 model
#Run first to export steady states to matrix file, then run
#reversible_steady_states.m in Matlab to generate plots
#Vary initial tBIM/tBid concentration

with(Statistics);
with(ArrayTools);
with(LinearAlgebra);
with(StringTools);

computeSteadyStates := proc(activator)
	global k1a, k1d, k2a, k2d, k3a, k3d, k4a, k4d, k5a, k6a, k6d, M0, B0, T0, g, a, K1D, K2A, K2D, K3A, K3D, K4A, K4D, K5A, K6A, K6D, b4, b_hat, m, d, b, t_hat, m_hat, f, d_ss, m_hat_ss, t_hat_ss, b_hat_ss, m_ss, b_ss, t_ss, b4_ss, biol, tss, tss_, ss;	
	#Kinetic parameters
	if (Compare(activator, 'tBIM')) then	
		k1a = 10^(-6.00);		 		#Bak activation (1/nM/s)
		k1d = 10^(-6);  				#Bak deactivation (1/s)
		k2a = 10^(-3.19);  				#Bak:Mcl-1 assoc (1/nM/s)
		k2d = 10^(-3.19);  				#Bak:Mcl-1 dissoc (1/s)
		k3a = 10^(-4.58);			    #tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-4.78);			    #tBim:Mcl-1 dissoc (1/s)	
		k4a = 10^(-1.63);				#Bak dimerisation (1/nM/s)
		k4d = 10^(-1.28);				#Bak dimerisation dissociation (1/s)
		k5a = 10^(-3.13);				#Bak auto-activation (1/nM/s)
		k5d = 0;
		k6a = 10^(-2.79);				#Bak multimer association rate (1/nM/s)
		k6d = 10^(-2.02);				#Bak multimer dissociation rate (1/s)
	elif (Compare(activator, 'tBid')) then
		#Values from fitted parameters (tBid case)
		k1a = 10^(-3.52);		 		#Bak activation (1/nM/s)
		k1d = 10^(-6);  				#Bak deactivation (1/s)
		k2a = 10^(-3.19);  				#Bak:Mcl-1 assoc (1/nM/s)
		k2d = 10^(-3.19);  				#Bak:Mcl-1 dissoc (1/s)
		k3a = 10^(-4.18);			    #tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.38);			    #tBid:Mcl-1 dissoc (1/s)
		k4a = 10^(-1.63);				#Bak dimerisation (1/nM/s)
		k4d = 10^(-1.28);				#Bak dimerisation dissociation (1/s)
		k5a = 10^(-3.13);				#Bak auto-activation (1/nM/s)
		k5d = 0;
		k6a = 10^(-2.79);				#Bak multimer association rate (1/nM/s)
		k6d = 10^(-2.02);				#Bak multimer dissociation rate (1/s)
	else
		return(FAIL);
	end if:	

	M0 := 20;
	B0 := 10;
	g := M0/B0;
	a := T0/M0;

	K1D := k1d/k1a/T0;
	K2A := k2a*B0/k1a/T0;
	K2D := k2d/k1a/T0;
	K3A := k3a*M0/k1a/T0;
	K3D := k3d/k1a/T0;
	K4A := 2*k4a*B0/k1a/T0;
	K4D := k4d/k1a/T0;
	K5A := k5a*B0/k1a/T0;
	K6A := k6a*B0/k1a/T0;
	K6D := k6d/k1a/T0;
	K7A := k7a*B0/k1a/T0/4;
	K8A := k8a*B0/k1a/T0;
	K8D := k8d/k1a/T0;
	K9A := k9a*X0/k1a/T0;
	K9D := k9d/k1a/T0;

	b_hat := K2D*m_hat/K2A/m;
	m     := K3D*(1-t)/K3A/t;
	d     := K4A*b_hat*b_hat/K4D;
	b     := K1D*b_hat/(t+K5A*b_hat);
	b4    := K6A*d*d/K6D;
	t_hat := 1-t;
	m_hat := 1-m-a+a*t;
	f := b_hat + b + d + 0*b4 + g*m_hat;
	t+th+thx-1, m+mh+alpha*th-1, b+bh+d+b4+lambda*mh+nu*xh-1, x+xh+kappa*thx-1, K2A*bh*m-K2D*mh, K3A*t*m-K3D*th, K6A*d*d-K6D*b4, K4A*bh*bh-K4D*d, b*t-K1D*bh+K5A*b*bh, K8A*x*bh-K8D*xh, K9A*t*x-K9D*thx);
	

	d_ss := Vector[row]([]);
	m_hat_ss := Vector[row]([]);
	t_hat_ss := Vector[row]([]);
	b_hat_ss := Vector[row]([]);
	m_ss := Vector[row]([]);
	b_ss := Vector[row]([]);
	t_ss := Vector[row]([]);
	b4_ss := Vector[row]([]);

	biol := (x,y) -> if (Im(eval(x, t=y)) <> 0 or eval(x,t=y)<0 or eval(x,t=y)>1) then false else true fi:

	for T0 in `$`(1..300) do
		tss := []:
		tss_ := [solve(f=1, t)]:
		for ss in tss_ do
			#For each ss check it is biologically meaningful
			#That is, every concentration for every species must be real and between 0 and 1
			if (biol(b4, ss) and biol(t, ss) and biol(b, ss) and biol(m, ss) and biol(d, ss) and biol(b_hat, ss) and biol(t_hat, ss) and biol(m_hat, ss)) then
				tss := [op(tss), ss]:
			end if:
		end do:
		b_ss := Concatenate(2, b_ss, Vector[column](map(x->eval(b,t=x), tss))):
		t_ss := Concatenate(2, t_ss, Vector[column](tss)):
		m_ss := Concatenate(2, m_ss, Vector[column](map(x->eval(m,t=x), tss))):
		d_ss := Concatenate(2, d_ss, Vector[column](map(x->eval(d,t=x), tss))):
		b_hat_ss := Concatenate(2, b_hat_ss, Vector[column](map(x->eval(b_hat,t=x), tss))):
		t_hat_ss := Concatenate(2, t_hat_ss, Vector[column](map(x->eval(t_hat,t=x), tss))):
		m_hat_ss := Concatenate(2, m_hat_ss, Vector[column](map(x->eval(m_hat,t=x), tss))):
		b4_ss := Concatenate(2, b4_ss, Vector[column](map(x->eval(b4,t=x),tss))):
	end do:
	

	#Plot them in Matlab so look consistent with other plots...
	A := Concatenate(1, b_ss, t_ss, m_ss, d_ss, b_hat_ss, t_hat_ss, m_hat_ss);
	ExportMatrix(cat(`./scripts/reversible_steady_states-`,activator,`.dat`), A, target=Matlab);
	#plotsetup(cps, plotoutput=cat(`./images/demos/rev_steady_states-`,activator,`-1.ps`), plotoptions=`border,width=8in,height=6in`);
	#print(plots[display](LineChart(convert(d_ss, listlist), view=[0..300, 0..1], gridlines=true, labels=['T0', 'd']),symbolsize=1,axes=boxed));
	#plotsetup(cps, plotoutput=cat(`./images/demos/rev_steady_states-`,activator,`-2.ps`), plotoptions=`border,width=8in,height=6in`);
	#print(plots[display](LineChart(convert(b_ss, listlist), view=[0..300, 0..1], gridlines=true, labels=['T0', 'b']),symbolsize=1,axes=boxed));
	#plotsetup(cps, plotoutput=cat(`./images/demos/rev_steady_states-`,activator,`-3.ps`), plotoptions=`border,width=8in,height=6in`);
	#print(plots[display](LineChart(convert(t_ss, listlist), view=[0..300, 0..1], gridlines=true, labels=['T0', 't']),symbolsize=1,axes=boxed));
	#plotsetup(cps, plotoutput=cat(`./images/demos/rev_steady_states-`,activator,`-4.ps`), plotoptions=`border,width=8in,height=6in`);
	#print(plots[display](LineChart(convert(m_ss, listlist), view=[0..300, 0..1], gridlines=true, labels=['T0', 'm']),symbolsize=1,axes=boxed));
	#plotsetup(cps, plotoutput=cat(`./images/demos/rev_steady_states-`,activator,`-5.ps`), plotoptions=`border,width=8in,height=6in`);
	#print(plots[display](LineChart(convert(b_hat_ss, listlist), view=[0..300, 0..1], gridlines=true, labels=['T0', 'b_hat']),symbolsize=1,axes=boxed));
	#plotsetup(cps, plotoutput=cat(`./images/demos/rev_steady_states-`,activator,`-6.ps`), plotoptions=`border,width=8in,height=6in`);
	#print(plots[display](LineChart(convert(m_hat_ss, listlist), view=[0..300, 0..1], gridlines=true, labels=['T0', 'm_hat']),symbolsize=1,axes=boxed));
	#plotsetup(cps, plotoutput=cat(`./images/demos/rev_steady_states-`,activator,`-7.ps`), plotoptions=`border,width=8in,height=6in`);
	#print(plots[display](LineChart(convert(t_hat_ss, listlist), view=[0..300, 0..1], gridlines=true, labels=['T0', 't_hat']),symbolsize=1,axes=boxed));
	#plotsetup(default);
endproc;

#computeSteadyStates('tBid');
computeSteadyStates('tBIM');	

