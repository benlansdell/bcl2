function [fits r2 bakmcl1r2 tbimmcl1r2] = fit_exp4_multi_diff_indirect(x0, filename_tbim, filename_tbid, plotfit)
	% fit_exp4_multi	Fit Experiment 4 quantitative observations to simulation results according to
	%			specified model, activator and a range of kinetic parameters. Uses trust-region reflective
	%			non-linear minimisation. Estimates only unknown kinetic rates k1a (Bak activation by BH3-only),
	%			k4a (Bak homodimerisation), k5a (Bak auto-activation), k1d (Bak* deactivation rate), 
	%			k4d (Bak* un-dimerisation rate). Here fit to very rough estimates of protein concentration
	%			taken from densitometry data and estimated quasi-equilibrium concentrations reached at 80 minute
	%			point in experiment. Uses diffusion-limited model. Fits the diffusion constant, 10^D. 
	%
	% Usage:
	%			[fits, r2, bakmcl1r2, tbidmcl1r2] = fit_exp4_multi_diff_indirect(x0, filename_tbim, filename_tbid, plotfit)
	%
	% Input:
	%			x0 = initial parameters for k4a, k5a, k6a, k4d, k6d, D
	%			filename_tbim = base of image filename for tbim simulation
	%			filename_tbid = base of image filename for tbid simluation
	%			plotfit (optional) = 1 for plot and write fitted parameters, 0 otherwise. Default = 1.
	%
	% Output:
	%			fits = vector of fitted values. In order: k4a, k5a, k6a, k4d, k6d, D
	%			r2 = r^2, coefficient of determination. Measure of the goodness-of-fit, value between 1 and 0
	%			bakmcl1r2 = r^2 value for fit to Bak:Mcl-1 observations only
	%			tbidcml1r2 = r^2 value for fit to tBid:Mcl-1 observations only
	%
	% Example(s):
	%			[fits, r2, bakmcl1r2, tbidmcl1r2] = fit_exp4_multi_diff([-4 -4 -4 -4 -4 -6], './images/test1', './images/test2');
	
	observations = observations_exp4_multi();
	model = 'reversible_exp4_BIAcore_indirect';
	method = 'trust-region-reflective';

	ub = [-1, -1, -1, -1, -1, -4];
	lb = [-6, -6, -6, -6, -6, -9];
	tbim_params = [1 2 3 4 5 6];
	tbid_params = [1 2 3 4 5 6];
	time_points = 60*[0 10 20 30 40 50 60 70 80] + 1;

	problem = createOptimProblem('lsqnonlin', 'x0', x0, 'objective', @(x) residuals(x, observations, model, tbim_params, tbid_params), 'ub', ub, 'lb',...
					 lb, 'options', optimset('Algorithm', method, 'TolX', 1e-4));
	%'nonlcon', @(x) multimer_con(x, model, tbim_params, tbid_params), 
	%, 'Display', 'iter-detailed'
	
	%Perform a 'local' search
	%tic
	%[fits, resnorm, exitflag, output] = lsqnonlin(problem);	
	%toc

	%Perform a 'global' search
	matlabpool 3;
	n_starts = 200;
	ms = MultiStart('TolX', 1e-4, 'MaxTime', 80000, 'Display', 'iter', 'UseParallel', 'always');
	[fits, resnorm, flagg, outptg, manyminsg] = run(ms, problem, n_starts);
	matlabpool close;
	
	f = run_simulation(model, 'tBid', T_p_exp4('tBid'), 0:(81*60), fits(tbid_params));
	[r2 bakmcl1r2 tbimmcl1r2] = compute_r2(observations, resnorm, f, time_points);
	if (nargin < 4) | (plotfit == 1)
		plot_fits('tBim', model, fits(tbim_params), observations(1:6,:), filename_tbim);
		plot_fits('tBid', model, fits(tbid_params), observations(7:13,:), filename_tbid);
	end
end

function resid = residuals(x, observations, model, tbim_params, tbid_params)
	resid = [solve_ode('tBim', model, x(tbim_params));
	         solve_ode('tBid', model, x(tbid_params))] - observations;
end

function ss = ss_resid(x, observations, model, tbim_params, tbid_params)
	resid = [solve_ode('tBim', model, x(tbim_params));
	         solve_ode('tBid', model, x(tbid_params))] - observations;
	ss = sum(sum(resid.^2));
end

function [f, h1, h2, h3, h4, h5] = run_simulation(model, activator, T_p, t, x)
	[k, ic, names] = load_model(activator, model);
	k(4,1) = 10^x(1);
	k(5,1) = 10^x(2);
	k(6,1) = 10^x(3);
	k(4,2) = 10^x(4);
	k(6,2) = 10^x(5);
	D = 10^x(6);
	if (nargout > 1)
		[f, h1, h2, h3, h4, h5] = bcl2model_diff(k, ic, D, names, t, T_p, activator);
	else
		f = bcl2model_diff(k, ic, D, names, t, T_p, activator);
	end
end

function f_scaled = solve_ode(activator, model, x)
	%Let it each equilibrium first...
	t = (-50*60):(81*60);
	time_points = (50*60)+60*[0 10 20 30 40 50 60 70 80] + 1;
	f = run_simulation(model, activator, T_p_exp4(activator), t, x);
	
	cytcc = f(9, time_points);
	cytcm = f(5, time_points);	
	bakmcl1 = f(8, time_points);
	tbimmcl1 = f(7, time_points);
	bak = f(1, time_points);
	abak = f(6, time_points);
	dimbak = f(4, time_points);
	mulbak = f(10, time_points);
	
	f_scaled = [cytcc; cytcm; bakmcl1; abak; dimbak; mulbak; tbimmcl1];
	
	if (strcmp(activator, 'tBim'))
		f_scaled = f_scaled(1:6,:);
	end
end

function plot_fits(activator, model, x, observations, filename)
	t = 0:(120*60);
	time_points = [0 10 20 30 40 50 60 70 80];
 	[f, h1, h2, h3, h4, h5] = run_simulation(model, activator, T_p_exp4(activator), t, x);
 	
	colororder = get(0, 'DefaultAxesColorOrder');
	figure(h1);
	xlim([0 120]);
	%Bak:Mcl-1
	line(time_points, observations(3,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(3,:));
	%Bak multimer
	%line(time_points, observations(5,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(5,:));
	saveplot(h1, [filename '-bak.eps']);

	figure(h2);
	xlim([0 120]);
	if (strcmp(activator, 'tBid'))
		line(time_points, observations(7,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(7,:));
	end
	saveplot(h2, [filename '-activator.eps']);

	figure(h3);
	xlim([0 120]);
	%Bak:Mcl-1
	line(time_points, observations(3,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(3,:));
	if (strcmp(activator, 'tBid'))
		%tBid:Mcl-1
		line(time_points, observations(7,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(7,:));
	end
	saveplot(h3, [filename '-mcl-1.eps']);

	figure(h4);
	xlim([0 120]);
	line(time_points, observations(1,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(9,:));
	line(time_points, observations(2,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(10,:));
	saveplot(h4, [filename '-cytc.eps']);

	figure(h5);
	xlim([0 120]);
	saveplot(h5, [filename '-bclxl.eps']);
end