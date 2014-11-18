function [conc, h1, h2, h3, h4] = bcl2model (k, ic, names, t, T_p, activator)
	% bcl2model			Solve ODE of specified Bcl2 interaction model. Use load_model() to generate
	%					model parameters. Includes all interactions possibly modelled -- direct activation,
	%					Bak multimerisation. These can be 'switched off' by the appropriate parameter set
	%					requested from load_model().
	%
	% Usage:
	%					[conc, h1, h2, h3, h4] = bcl2model(k, ic, names, t, T_p, activator)
	%
	% Input:
	%					k = 8x2 matrix of kinetic parameters containing association and dissocation rates
	%						for each of the 8 reactions modelled
	%					ic = 10x1 vector of initial conditions for each of the 10 species modelled.
	%						In order of species listed in 'names' cell provided by load_model()
	%					names = 10x1 cell array of names of species. Provided by load_model()
	%					t = vector of time points to calculate solution for
	%					T_p = BH3-only spike-in function. Provided by T_P_*() methods
	%					activator = one of 'tBim' or 'tBid'. BH3-only stimulus protein
	%
	% Output:
	%					conc = 10xs matrix of concentrations, where s denotes the length of input vector 't'.
	%							Concentration of each species at each timepoint.
	%					h1 = handle for plot object for Mcl-1 and Mcl-1 complexes plot
	%					h2 = handle for plot object for Bak and Bak complexes plot
	%					h3 = handle for plot object for tBid/tBim and tBid/tBim complexes plot
	%					h4 = handle for plot object for cytochrome c complexes
	%
	% Examples:
	%					[k, ic, names] = load_model('tBid', 'reversible_exp4');
	%					T_p = T_p_exp4('tBid');
	%					[conc, h1, h2, h3, h4] = bcl2model(k, ic, names, 0:5400, T_p, 'tBid');
	%					saveplot(h1, './images/test.eps');
	warning('bcl2model_comp depracated. Do not use');

	%'Cooperation coefficient' for Bak multimerisation
	n = 2;
	%Reactions
	%Activation of Bak
	v1 = @(x) k(1,1)*x(1)*x(2) - k(1,2)*x(6);
	%Formation of Bak:Mcl-1 complex
	v2 = @(x) k(2,1)*x(3)*x(6) -k(2,2)*x(8);
	%Formation of tBim(d)/Mcl-1 complex
	v3 = @(x) k(3,1)*x(2)*x(3) - k(3,2)*x(7);
	%Dimerisation of Bak*
	v4 = @(x) k(4,1)*x(6)*x(6) - k(4,2)*x(4);
	%Auto-activation of Bak
	v5 = @(x) k(5,1)*x(6)*x(1);
	%Multimerisation of Bak*:Bak*
	v6 = @(x) k(6,1)*x(4).^n - k(6,2)*x(10);
	%Cytocrhome C release
	v7 = @(x) k(7,1)*x(10)*x(5);
	%Compartment model additions:
	%Transport of Mcl-1
	v8 = @(x) k(8,1)*x(11) - k(8,2)*x(3);
	%Formation of tBim(d):Mcl-1_c complex
	v9 = @(x) k(9,1)*x(11)*x(12)-k(9,2)*x(13);
	%Transport of tBim(d)_c
	v10 = @(x) k(8,1)*x(12) - k(8,1)*x(2);
	%Transport of tBim(d):Mcl-1_c
	v11 = @(x) k(8,1)*x(13) - k(8,2)*x(7);

	%Derivative function for solving IVP problem
	df = @(t, x) [	-v1(x) - v5(x);			v10(x) - v3(x) + T_p(t); 
					v8(x) - v2(x) - v3(x); 	v4(x) - n*v6(x);
					-v7(x);					v1(x) + v5(x) - v2(x) - 2*v4(x);
					v11(x) + v3(x);			v2(x);
					v7(x);					v6(x);
					-v8(x) - v9(x);			-v10(x) - v9(x);
					v9(x) - v11(x)];

	%Needs to be set or may skip over tBid/tBim being added...
	options = odeset('MaxStep', 5, 'NonNegative', [1 1 1 1 1 1 1 1 1 1 1 1 1]);
	solver = @ode15s;
	[T, f] = solver(df, t, ic, options);
	try
		assert(length(T) == length(t))
	catch ME
		warning(['ODE solver failed with parameters: ' num2str(reshape(k, 1, prod(size(k))))]);
		T = t;
		f = repmat([0], length(t), length(ic));
	end

	%Return result and plots
	f = f';
	bakmcl1 = f(8,:);
	abak = f(6,:);
	bak = f(1,:);
	mcl1 = f(3,:);
	tbimmcl1 = f(7,:);
	tbim = f(2,:);
	dimbak = f(4,:);
	cytcC = f(9,:);
	cytcM = f(5,:);
	mulbak = f(10,:);

	colororder = get(0, 'DefaultAxesColorOrder');
	h1 = figure();
	plot(t/60, bak, 'Color', colororder(1,:));
	line(t/60, abak, 'Color', colororder(2,:));
	line(t/60, bakmcl1, 'Color', colororder(3,:));
	line(t/60, dimbak, 'Color', colororder(4,:));
	line(t/60, mulbak, 'Color', colororder(5,:));
	legend('Bak', 'Bak*', 'Bak:Mcl-1', 'Dim. Bak*', 'Multi. Bak*', 'Location', 'NorthEast');
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 12]);
	xlim([0 length(t)/60]);
	
	h2 = figure();
	plot(t/60, tbim, 'Color', colororder(6,:));
	line(t/60, tbimmcl1, 'Color', colororder(7,:));
	if (strcmp(activator, 'tBid'))
		legend(activator, [activator ':Mcl-1'], 'Location', 'NorthEast');
	else
		legend(activator, [activator ':Mcl-1'], 'Location', 'SouthEast');
	end
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 23]);
	xlim([0 length(t)/60]);
	
	h3 = figure();
	plot(t/60, mcl1, 'Color', colororder(8,:));
	line(t/60, bakmcl1, 'Color', colororder(3,:))
	line(t/60, tbimmcl1, 'Color', colororder(7,:))
	if (strcmp(activator, 'tBid'))
		legend('Mcl-1', 'Bak*:Mcl-1', [activator ':Mcl-1'], 'Location', 'NorthEast');
	else
		legend('Mcl-1', 'Bak*:Mcl-1', [activator ':Mcl-1'], 'Location', 'SouthEast');
	end
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 23]);
	xlim([0 length(t)/60]);
	
	h4 = figure();
	plot(t/60, cytcC, 'Color', colororder(9,:));
	line(t/60, cytcM, 'Color', colororder(10,:))
	legend('Cyt c S/N', 'Cyt c mito', 'Location', 'NorthEast');
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 6]);
	xlim([0 length(t)/60]);
	
	%Return concentrations of experiment
	conc = f;
end
