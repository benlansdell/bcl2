function R = robust_analysis_bentele_6(model, activator, filename)
	% robust_analysis_bentele		Function to analyse the sensitivity of model defined by specified
	%						kinetic parameters. Use 'global' method from Bentele 04. 
	%
	% Usage:
	%
	%						R = robust_analysis_bentele_6(model, activator, filename)
	%
	% Input:
	%						k = 9x2 matrix of kinetic parameters
	%						ic = 13x1 vector of initial concentrations
	%						directact = 1 (default) for direct activation and 0 for indirect
	%						filename = output EPS file for 2d bar plot
	%
	% Output:
	%						R = structure with the following components:
	%						R.S = Sensitivity matrix
	% 
	% Examples:
	%						R = robust_analysis_bentele_6('reversible_exp14_both', 'tBid', './images/test.eps');

	dk = 0.1;

	%Direct activation model
	x = [1 1 2 2 3 3 4 4 5 6  6 7];
	y = [1 2 1 2 1 2 1 2 1 1  2 1];
	z = [1 2 3 4 5 6 7 8 9 10 11 12 13];
	[k, ic, names] = load_model14(activator, model);
	t = (-50*60):(120*60);
	f = bcl2model(k, ic, names, t, T_p_exp14(activator), activator);
	S = zeros(13,12);
	for i=1:13
		for j=1:12
			i
			j
			k_ = k;
			k_(x(j), y(j)) = k(x(j), y(j))*dk;
			f_ = bcl2model(k_, ic, names, t, T_p_exp14(activator), activator);
			c = f(z(i),:);
			c_ = f_(z(i),:);
			C = sum(c, 2);
			C_ = sum(abs(c_-c),2);
			S(i,j) = log10(k(x(j),y(j)))*C_/C/log10(dk);
		end
	end
	R = struct('S', S);
	
	colororder = get(0, 'DefaultAxesColorOrder');
	bar3(S, 0.15);
	ylabel('Species');
	xlabel('Parameter');
	set(gca,'YTickLabel', names)
	set(gca,'XTickLabel', {'k1a', 'k1d', 'k2a', 'k2d', 'k3a', 'k3d', 'k4a', 'k4d', 'k5a', 'k6a', 'k6d', 'k7a'})
	set(gca,'Zlim',[0 16]);
	%hmo = HeatMap(S, 'RowLabels', names, 'ColumnLabels', {'k1a', 'k1d', 'k2a', 'k2d', 'k3a', 'k3d', 'k4a', 'k4d', 'k5a', 'k6a', 'k6d', 'k7a'});
	%h = addXLabel(hmo, 'Parameters');
	%g = addYLabel(hmo, 'Species');
	%plot(hmo);
	set(gca, 'FontSize', 6);
	saveplot(gcf, filename);
end
