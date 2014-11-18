function R = robust_analysis(k, ic, observations)
	% robust_analysis		Function to analyse the 'robustness' of model defined by specified
	%						kinetic parameters. At steady-state compute community matrix, relative
	%						change in concentration for change in kinetic parameter, change in 
	%						quality of least-squares fit for change in kinetic parameter -- discovery
	%						of 'sloppy' parameters.
	%
	% Usage:
	%
	%						R = robust_analysis(k, ic, observations)
	%
	% Input:
	%						k = 7x2 matrix of kinetic parameters
	%						ic = 10x1 vector of initial concentrations
	%						observations = 13x9 matrix of 'densitometry' which was used to fit model to
	%
	% Output:
	%						R = structure with the following components:
	%						R.is_stable = 1 if steady-state is stable, 0 if unstable
	%						R.ss = 10x1 vector of steady-state concentrations
	%						R.A = community matrix about steady state
	%						R.dX = Jacobian of steady-state concentration function X(k)
	%
	%						R.dR2 = Gradient of R^2 (coefficient of determination)
	% 
	% Examples:
	%						[k, ic] = load_model('tBid', 'reversible_exp4_BIAcore');
	%%						obs = observations_exp4_multi();
	%						R = robust_analysis(k, ic, obs);

	dk = 0.1;
	x = [1 1 2 2 3 3 4 4 5 6 6 7];
	y = [1 2 1 2 1 2 1 2 1 1 2 1];
	z = [1 2 3 4 5 6];
	ss = steadystate_conc(k, ic);
	[is_stable, A] = compute_stability(k, ic);
	dX = zeros(6,12);
	for i=1:6
		for j=1:12
			i
			j
			k_ = k;
			k_(x(j), y(j)) = k(x(j), y(j))+dk;
			ss_ = steadystate_conc(k_, ic);
			dX(i,j) = (ss_(z(i)) - ss(z(i)))/dk/ss(z(i)); 
		end
	end
	[U,S,V] = svd(dX);
	R = struct('ss', ss, 'is_stable', is_stable, 'A', A, 'dX', dX, 'U', U, 'S', S, 'V', V);
end
