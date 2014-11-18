function R = robust_analysis_direct(k, ic, observations)
	% robust_analysis		Function to analyse the 'robustness' of model defined by specified
	%						kinetic parameters. At steady-state compute community matrix, relative
	%						change in concentration for change in kinetic parameter, change in 
	%						quality of least-squares fit for change in kinetic parameter -- discovery
	%						of 'sloppy' parameters.
	%
	% Usage:
	%
	%						R = robust_analysis_direct(k, ic, observations)
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
	%						R = robust_analysis_direct(k, ic, obs);

	dk = 0.1;

	x = [1 1 2 2 3 3 4 4 5 6  6];
	y = [1 2 1 2 1 2 1 2 1 1  2];
	z = [1 2 3 4 5 6];
	a = [1 1 3 3 2 2 6 6 6 4  4];
	b = [2 2 6 6 3 3 6 6 1 4  4];
	c = [6 6 8 8 7 7 4 4 1 10 10];
	ss = steadystate_conc(k, ic, p, u);
	[is_stable, A] = compute_stability(k, ic);
	dX = zeros(6,11);
	for i=1:6
		for j=1:11
			i
			j
			k_ = k;
			k_(x(j), y(j)) = k(x(j), y(j))+dk;
			ss_ = steadystate_conc(k_, ic);
			v = k(x(j),1)*ss(a(j))*ss(b(j)) - k(x(j),2)*ss(c(j))
			v_ = k_(x(j),1)*ss_(a(j))*ss_(b(j)) - k_(x(j),2)*ss_(c(j))
			dv = v_-v;
			dX(i,j) = v*(ss_(z(i)) - ss(z(i)))/dv/ss(z(i)); 
		end
	end
	[U,S,V] = svd(dX);
	R = struct('ss', ss, 'is_stable', is_stable, 'A', A, 'dX', dX, 'U', U, 'S', S, 'V', V);
end
