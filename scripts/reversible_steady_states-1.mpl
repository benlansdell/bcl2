#Script to compute the steady state(s?) of the reversible Bcl2 model
#Run first to export steady states to matrix file, then run
#reversible_steady_states.m in Matlab to generate plots
#Vary initial tBIM/tBid concentration

with(Statistics);
with(ArrayTools);
with(LinearAlgebra);
with(StringTools);
with(VectorCalculus);

computeStability := proc(bss, tss, mss, dss, b_hatss, t_hatss, m_hatss, gss)
	t_hats := 1-ts;
	m_hats := 1-ms-a+a*ts;
	b_hats := 1-bs-ds-l*m_hats;
	S := Matrix([[-1,0,0,0,-1,0],[0,0,-1,0,0,0],[0,-1,-a,0,0,0],[0,0,0,1,0,0],[0,0,0,0,0,-1]]);
	J := Vector([bs*ts - K1D*b_hats,K2A*ms*b_hats-K2D*m_hats,K3A*ts*ms-K3D*t_hats,K4A*b_hats*b_hats-K4D*ds,K5A*bs*b_hats,K6A*ds*gs]);
	j := Jacobian(J, [bs,ts,ms,ds,gs]=[bss,tss,mss,dss,gss]);
	A := Multiply(S,j);
	evals := Eigenvalues(A);
	real_negative := (x) -> if (Im(x) = 0 and Re(x) < 0) then true else false fi;
	retval_b := real_negative(evals[1]) and real_negative(evals[2]) and real_negative(evals[3]) and real_negative(evals[4]) and real_negative(evals[5]);
	retval_i := (x) -> if (x = true) then 1 else 0 fi;
	return(retval_i(retval_b));
end proc;

computeSteadyStates := proc(activator)
	global k1a, k1d, k2a, k2d, k3a, k3d, k4a, k4d, k5a, M0, B0, T0, l, a, K1D, K2A, K2D, K3A, K3D, K4A, K4D, K5A, K6A, b_hat, m, d, b, t_hat, m_hat, f, d_ss, m_hat_ss, t_hat_ss, b_hat_ss, m_ss, b_ss, t_ss, biol, tss, tss_, ss;	
	if (Compare(activator, 'tBIM')) then	
		#Values from fitted parameters (tBIM case)
		k1a := 10^(-6);
		k1d := 10^(-6);
		k2a := 10^(-2.55);
		k2d := 10^(-2.77);
		k3a := 10^(-4.7);
		k3d := 10^(-4.94);
		k4a := 10^(-1.14);
		k4d := 10^(-1.52);
		k5a := 10^(-3.13);
		k6a := 10^(-4);
	elif (Compare(activator, 'tBid')) then
		#Values from fitted parameters (tBid case)
		k1a := 10^(-3.7800);
		k3a := 10^(-4.2966);
		k3d := 10^(-2.5351);
		k2a := 10^(-2.5500);
		k4a := 10^(-1.1400);
		k5a := 10^(-3.1300);
		k1d := 10^(-6.0000);
		k2d := 10^(-2.7700);
		k4d := 10^(-1.5200);
		k6a := 10^(-4);
	else
		return(FAIL);
	end if:	
	
	M0 := 20;
	B0 := 10;
	l := M0/B0;
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

	b_hat := K2D*m_hat/K2A/m;
	m     := K3D*(1-t)/K3A/t;
	d     := K4A*b_hat*b_hat/K4D;
	b     := K1D*b_hat/(t+K5A*b_hat);
	t_hat := 1-t;
	m_hat := 1-m-a+a*t;
	f := b_hat + b + d + l*m_hat;

	d_ss := Vector[row]([]);
	m_hat_ss := Vector[row]([]);
	t_hat_ss := Vector[row]([]);
	b_hat_ss := Vector[row]([]);
	m_ss := Vector[row]([]);
	b_ss := Vector[row]([]);
	t_ss := Vector[row]([]);
	stab_ss := Vector[row]([]);

	biol := (x,y) -> if (Im(eval(x, t=y)) <> 0 or eval(x,t=y)<0 or eval(x,t=y)>1) then false else true fi;

	for T0 in `$`(1..300) do
		tss := []:
		tss_ := [solve(f=1, t)]:
		for ss in tss_ do
			#For each ss check it is biologically meaningful
			#That is, every concentration for every species must be real and between 0 and 1
			if (biol(t, ss) and biol(b, ss) and biol(m, ss) and biol(d, ss) and biol(b_hat, ss) and biol(t_hat, ss) and biol(m_hat, ss)) then
				tss := [op(tss), ss]:
			end if:
		end do:
		#Don't proceed with more than one valid steady-state...
		if (nops(tss) <> 1) then return(FAIL) end if;
		b_ss := Concatenate(2, b_ss, Vector[column](map(x->eval(b,t=x), tss)));
		t_ss := Concatenate(2, t_ss, Vector[column](tss)):
		m_ss := Concatenate(2, m_ss, Vector[column](map(x->eval(m,t=x), tss))):
		d_ss := Concatenate(2, d_ss, Vector[column](map(x->eval(d,t=x), tss))):
		b_hat_ss := Concatenate(2, b_hat_ss, Vector[column](map(x->eval(b_hat,t=x), tss))):
		t_hat_ss := Concatenate(2, t_hat_ss, Vector[column](map(x->eval(t_hat,t=x), tss))):
		m_hat_ss := Concatenate(2, m_hat_ss, Vector[column](map(x->eval(m_hat,t=x), tss))):
		stab_ss := Concatenate(2, stab_ss, Vector[column]([computeStability(eval(b, t=tss[1]), tss[1], eval(m, t=tss[1]), eval(d, t=tss[1]), eval(b_hat, t=tss[1]), eval(t_hat, t=tss[1]), eval(m_hat, t=tss[1]), 1)])):
	end do;
	
	#Plot them in Matlab so look consistent with other plots...
	A := Concatenate(1, b_ss, t_ss, m_ss, d_ss, b_hat_ss, t_hat_ss, m_hat_ss, stab_ss);
	ExportMatrix(cat(`./scripts/reversible_steady_states-`,activator,`-1.dat`), Matrix(A), target=Matlab);
end proc;

#computeSteadyStates('tBid');
computeSteadyStates('tBIM');	

