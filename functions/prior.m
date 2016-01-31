function lp = prior(k)
	%Evaluate prior probability given a set of proposed parameters
	%The prior considers knowledge obtained about the kinetics
	%obtained before the experiments were performed -- i.e. data from BIAcore
	%measurements.

	BIAcore_a3 = log10(8.10/3.20);
	BIAcore_d3 = log10(752/2.96);
	BIAcore_a9 = log10(1.55/3.86);
	BIAcore_d9 = log10(39.2/476);

	%Quantified intial conditions for Bak and Bcl-xl
	mu_bak = 0.89;
	mu_bclxl = 1.15;
	sigma_bak = 0.13;
	sigma_bclxl = 0.11;

	ub = 1;
	lb = -8;

	z_logrates = [k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8), k(9), k(10), k(11), k(12), k(13), ...
	k(14), k(15), k(16), k(17), k(18), k(19), k(20), k(21), k(22), k(23)];
	z_conc = [k(24), k(25)];

	if any(z_logrates < lb | z_logrates > ub)
		lp = -Inf;
		return;
	end

	if any(z_conc <= 0)
		lp = -Inf;
		return;
	end

	mu_zero = -7;
	mu = -3.5;
	%Flat prior for parameters without BIAcore knowledge
	s1 = 2;
	%Moderate prior for parameters with BIAcore knowledge
	s2 = 0.5;
	%Sharp prior for parameters we wish to set to effectively zero
	s3 = 0.1;

	%List of parameters
	%k1a_tbim, k1a_tbid, k2a, k3a_tbim, k3a_tbid, k4a, k5a, k6a, k7a, k8a, k9a_tbim, k9a_tbid, k1d, k2d,
	%k3d_tbim, k3d_tbid, k4d, k5d, k6d, k7d, k8d, k9d_tbim, k9d_tbid, IC_mu_bak, IC_mu_bclxl

	%Things being set to (effective) zero: k5d, k7d

		%k1a_tbim, k1a_tbid, k2a,  k3abim,k3abid, k4a,  k5a,  k6a,  k7a,  k8a,   k9abim,k9abid, k1d, 
	%k2d,  k3dbim,k3dbid, k4d,   k5d,   k6d,   k7d,   k8d,   k9dbim,k9dbid, k3abidbim  k9abidbim
	%k3dbidbim   k9dbidbim
	z = [k(1),     k(2),     k(3), k(4),  k(5),   k(6), k(7), k(8), k(9), k(10), k(11), k(12),  k(13), ...
	k(14), k(15), k(16),  k(17), k(18), k(19), k(20), k(21), k(22), k(23),  k(5)-k(4), k(12)-k(11), ...
	k(16)-k(15), k(23)-k(22), k(24), k(25)];

		%   k1a_tbim, k1a_tbid, k2a,  k3abim,k3abid, k4a,  k5a,  k6a,  k7a,  k8a, k9abim,k9abid, k1d, 
	%k2d,   k3dbim,k3dbid, k4d,   k5d,   k6d,   k7d,   k8d, k9dbim, k9dbid, k3abidbim,  k9abidbim, k3dbidbim,  k9dbidbim,   IC_bak, IC_bclxl
	z_mu = [mu,       mu,       mu,   mu,    mu,     mu,   mu,   mu,   mu,   mu,  mu,    mu,     mu,  ...
	 mu,    mu,    mu,     mu,  mu_zero, mu,  mu_zero, mu,  mu,     mu,     BIAcore_a3, BIAcore_d3, BIAcore_a9, BIAcore_d9, mu_bak, mu_bclxl];

	   %      k1a_tbim, k1a_tbid, k2a,  k3abim,k3abid, k4a,  k5a,  k6a,  k7a,  k8a,   k9abim, k9abid, k1d, 
	%k2d,   k3dbim, k3dbid, k4d,   k5d,   k6d,   k7d,   k8d, k9dbim, k9dbid, k3abidbim  k9abidbim, k3dbidbim  k9dbidbim, IC_bak,    IC_bclxl]
	stddev = [s1,       s1,       s1,   s1,    s1,     s1,   s1,   s1,    s1,  s1,    s1,     s1,     s1, ...
	 s1,    s1,     s1,     s1,    s3,    s1,    s3,    s1,  s1,     s1,     s2         s2         s2         s2,        sigma_bak, sigma_bclxl];

	p = normpdf(z, z_mu, stddev);
	lp = sum(log(p));
end