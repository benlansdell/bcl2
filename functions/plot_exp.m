function plot_exp(model, activator, filename)
	% plot_exp		Plot concentration vs. time for specified Experiment 4 model 
	%				and kinetic parameters.
	%
	% Usage:
	%			plot_exp(func, model, activator, filename)
	%
	% Input:
	%			model = determines kinetic parameters, initial conditions and other parameters.
	%					See help load_model for options.
	%			activator = one of 'tBid' or 'tBim'
	%			filename = base for filenames to save plots. Plots are saved to
	%						'filename'-mcl1.eps, 'filename'-bak.eps, 'filename'-activator.eps
	%
	% Example(s):
	%			plot_exp('reversible_exp4', 'tBid', './images/test');
	
	[k, ic, names, D] = load_model(activator, model);
	T_p = T_p_exp4(activator);
	T = 60*120;
	if (D ~= 0)
		[conc, h1, h2, h3, h4] = bcl2model_diff(k, ic, names, 0:T, T_p, activator);
	else
		[conc, h1, h2, h3, h4] = bcl2model(k, ic, names, 0:T, T_p, activator);
	end
	saveplot(h1, [filename '-bak.eps']);
	saveplot(h2, [filename '-activator.eps']);
	saveplot(h3, [filename '-mcl1.eps']);
	saveplot(h4, [filename '-cytc.eps']);
end