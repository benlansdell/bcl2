function R = robust_analysis_bentele_rndparam_6(n, basename)
	% robust_analysis_bentele_rndparam_6		Function to analyse the sensitivity of model defined by specified
	%						kinetic parameters. Picks parameters randomly and collects statistics
	%						on sensitivity.  Uses 'global' method from Bentele 04. 
	%
	% Usage:
	%
	%						R = robust_analysis_bentele_rndparam_6(n, basename)
	%
	% Input:
	%						n = number of random parameters to generate
	%						basename = output EPS base filename for 2d bar plots of mean and std deviation
	%
	% Output:
	%						R = structure with the following components:
	%						R.mu = matrix of mean sensitivities
	%						R.sigma = matrix of standard deviation for sensitivities
	% 
	% Examples:
	%						R = robust_analysis_bentele_rndparam_6(1000, './images/test');

	%Change parameter by 10%
	dk = 0.1;
	activator = 'tBid';
	model = 'random_parameters_both';

	%Direct activation model
	x = [1 1 2 2 3 3 4 4 5 6  6 7];
	y = [1 2 1 2 1 2 1 2 1 1  2 1];
	z = [1 2 3 4 5 6 7 8 9 10 11 12 13];
	[k, ic, names] = load_model6(activator, model);
	t = (-50*60):(120*60);
	f = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
	S = zeros(n,13,12);
	
	matlabpool 4;
	parfor m = 1:n
		m
		[k, ic, names] = load_model6(activator, model);
		f = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
		S_ = zeros(1,13,12);
		for i=1:13
			for j=1:12
				k_ = k;
				k_(x(j), y(j)) = k(x(j), y(j))*dk;
				f_ = bcl2model(k_, ic, names, t, T_p_exp6(activator), activator);
				c = f(z(i),:);
				c_ = f_(z(i),:);
				C = sum(c, 2);
				C_ = sum(abs(c_-c),2);
				S_(1,i,j) = -log10(k(x(j),y(j)))*C_/C/dk;
			end
		end
		S(m,:,:) = S_(1,:,:);
	end
	matlabpool close;
	R = struct('mu', reshape(mean(S,1),13,12), 'sigma', reshape(std(S,0,1),13,12));
	colororder = get(0, 'DefaultAxesColorOrder');
	bar3(R.mu, 0.15);
	ylabel('Species');
	xlabel('Parameter');
	set(gca,'YTickLabel', names)
	set(gca,'XTickLabel', {'k1a', 'k1d', 'k2a', 'k2d', 'k3a', 'k3d', 'k4a', 'k4d', 'k5a', 'k6a', 'k6d', 'k7a'});
	set(gca,'Zlim', [0 200]);
	set(gca, 'FontSize', 6);
	saveplot(gcf, [basename '_mu.eps']);

        bar3(R.sigma, 0.15);
        ylabel('Species');
        xlabel('Parameter');
        set(gca,'YTickLabel', names)
        set(gca,'XTickLabel', {'k1a', 'k1d', 'k2a', 'k2d', 'k3a', 'k3d', 'k4a', 'k4d', 'k5a', 'k6a', 'k6d', 'k7a'});
	set(gca,'Zlim', [0 200]);
        set(gca, 'FontSize', 6);
        saveplot(gcf, [basename '_sigma.eps']);	
      
	save('./scripts/robust_analysis_bentele_rndparam_6_Smatrix.mat', 'S')
end
