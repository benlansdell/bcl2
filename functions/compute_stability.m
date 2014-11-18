function [stability, community, evals] = compute_stability(k, ic)
	%compute_stability		Function to compute the stability of a steady state computed using
	%						provided kinetic parameters. The linearised 'community' matrix is 
	%						also returned. Reduces 
	%
	% Usage:
	%
	% Input:
	%
	% Output:
	%
	% Examples:
	%						[k, ic] = load_model('tBid', 'reversible_exp4');
	%						compute_stability(k, ic);

	B0 = ic(1);
	T0 = ic(2);
	M0 = ic(3);
	C0 = ic(5);

	a = T0/M0;
	l = M0/B0;
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

	syms b t m d g bh th mh gh b4;
	%Stoichiometric matrix
	S = [[-1, 0, 0, 0,-1, 0, 0];
		 [ 0, 0,-1, 0, 0, 0, 0];
		 [ 0,-1,-a, 0, 0, 0, 0];
		 [ 0, 0, 0, 1, 0,-1, 0];
		 [ 0, 0, 0, 0, 0, 0,-1];
		 [ 1,-l, 0,-1, 1, 0, 0]];
	%Rate vector
	J = [b*t - K1D*bh, K2A*m*bh-K2D*mh, K3A*t*m-K3D*th, K4A*bh*bh-K4D*d, K5A*b*bh, K6A*d*d-K6D*b4, K7A*b4*g];
	%Substitute conservation relations
	J = subs(J, {mh, th, b4}, {1-m-a*(1-t),1-t,1-b-bh-d-l*(1-m-a*(1-t))});
	%Compute jacobian
	j = jacobian(J, [b, t, m, d, g, bh]);
	%Evaluate at steady state
	[bs, ts, ms, ds, gs, bhs] = steadystate_conc(k, ic);
	j = subs(j, {b, t, m, d, g, bh}, {bs, ts, ms, ds, gs, bhs});
	%Form linearised system about steady-state
	community = S*j;
	%Compute eigenvalues
	evals = eig(community);
	stability = all(evals<0);
end