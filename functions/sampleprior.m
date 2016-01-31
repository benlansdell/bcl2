function rd = sampleprior()
	%Evaluate prior probability given a set of proposed parameters
	%The prior considers knowledge obtained about the kinetics
	%obtained before the experiments were performed -- i.e. data from BIAcore
	%measurements.

	%warning('Not implemented yet')

	BIAcore_a3 = log10(8.10/3.20);
	BIAcore_d3 = log10(752/2.96);
	BIAcore_a9 = log10(1.55/3.86);
	BIAcore_d9 = log10(39.2/476);


	%Quantified intial conditions for Bak and Bcl-xl
	mu_bak = 0.89;
	mu_bclxl = 1.15;
	sigma_bak = 0.13;
	sigma_bclxl = 0.11;

	mu_zero = -7;
	mu = -3.5;
	%Flat prior for parameters without BIAcore knowledge
	s1 = 2;
	%Moderate prior for parameters with BIAcore knowledge
	s2 = 0.5;
	%Sharp prior for parameters we wish to set to effectively zero
	s3 = 0.1;

	ub = 1;
	lb = -8;

	%List of parameters
	%k1a_tbim, k1a_tbid, k2a, k3a_tbim, k3a_tbid, k4a, k5a, k6a, k7a, k8a, k9a_tbim, k9a_tbid, k1d, k2d,
	%k3d_tbim, k3d_tbid, k4d, k5d, k6d, k7d, k8d, k9d_tbim, k9d_tbid

	%Things being set to (effective) zero: k5d, k7d

		%k1a_tbim, k1a_tbid, k2a,  k4a,  k5a,  k6a,  k7a,  k8a,   k1d, 
	%k2d,  k4d,   k5d,   k6d,   k7d,   k8d,   k3abidbim  k9abidbim
	%k3dbidbim   k9dbidbim
	%z = [k(1),     k(2),     k(3), k(6), k(7), k(8), k(9), k(10), k(13), ...
	%k(14), k(17), k(18), k(19), k(20), k(21), k(5)-k(4), k(12)-k(11), ...
	%k(16)-k(15), k(23)-k(22)];

		%   k1a_tbim, k1a_tbid, k2a,  k4a,  k5a,  k6a,  k7a,  k8a,   k1d, 
	%k2d,   k4d,   k5d,   k6d,   k7d,   k8d, k3abidbim,  k9abidbim, k3dbidbim,  k9dbidbim
	z_mu = [mu,       mu,       mu,   mu,   mu,   mu,   mu,   mu,    mu,...
	 mu,    mu,  mu_zero, mu,  mu_zero, mu, BIAcore_a3, BIAcore_d3, BIAcore_a9, BIAcore_d9, mu_bak, mu_bclxl];

	   %      k1a_tbim, k1a_tbid, k2a,  k4a,  k5a,  k6a,  k7a,  k8a,   k1d, 
	%k2d,   k4d,   k5d,   k6d,   k7d,   k8d, k3abidbim  k9abidbim, k3dbidbim  k9dbidbim]
	stddev = [s1,       s1,       s1,   s1,   s1,   s1,   s1,   s1,    s1,...
	 s1,    s1,    s3,    s1,    s3,    s1,  s2         s2         s2         s2,        sigma_bak, sigma_bclxl];

	validvalues = 0;
	while (validvalues == 0)
		validvalues = 1;
		%Sample from p
		nz = size(z_mu, 2);
		zs = randn(1, nz);
		s = z_mu + zs.*stddev;
		zs = randn(1, 4);
		t = z_mu(1:4) + zs.*stddev(1:4);
	
		%
		%[s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), s(9), ...
		%                                          %%     %%
		%s(10), s(11), s(12), s(13), s(14), s(15), s(16), s(17), ...
		%%      %%
		%s(18), s(19)];
	
		rd = [s(1), s(2), s(3), t(1)+s(16), t(1), s(4), s(5), s(6), s(7), s(8), t(2)+s(17), ...
		t(2), s(9), s(10), t(3)+s(18), t(3), s(11), s(12), s(13), s(14), s(15), t(4)+s(19), t(4), s(20), s(21)];
	
		z_logrates = [rd(1), rd(2), rd(3), rd(4), rd(5), rd(6), rd(7), rd(8), rd(9), rd(10), rd(11), rd(12), rd(13), ...
		rd(14), rd(15), rd(16), rd(17), rd(18), rd(19), rd(20), rd(21), rd(22), rd(23)];
		z_conc = [rd(24), rd(25)];
	
		if any(z_logrates < lb | z_logrates > ub)
			validvalues = 0;
		end
		if any(z_conc <= 0)
			validvalues = 0;
		end
	end 
end