function optimisationtest()
	%optimisationtest		Script to test speed and accuracy of different optimisation methods
	%						in solving non-linear least squares problem of parameter estimation...
	%						Fit a number of models with different number of free parameters to fit...
	%						In each case only use one set of initial parameters. Use 'typical' kinetic
	%						values: 10^(-4)
	%
	%Output:
	%
	%Initial conditions: -4 -4 -4...
	%
	%Parameters:
	%--------------------
	%Base fit obtained: 			-4.4376      -3.0841	-3.5274		-3.4812 																					 with norm 0.44999 in 209.9964 seconds using trust-region-reflective
	%Base fit obtained:			-3.1578      -3.05769	-3.62706	-12.0644 																					 with norm 0.4602 in 161.8628 seconds using levenberg-marquardt
	%Diffusion fit obtained: 	-4.4291      -3.1457	-3.5126		-3.4759		-4.9492																			 with norm 0.44927 in 230.9859 seconds using trust-region-reflective
	%Diffusion fit obtained: 	-3.08083     -2.89032	-3.66807	-14.1813	-7.6758 																		 with norm 0.46616 in 169.4869 seconds using levenberg-marquardt
	%All-forwards fit obtained: 	-4.7573      -3.526		-3.229		-3.8645		-3.1378		-3.1268		-3.7598													 with norm 0.72759 in 215.1207 seconds using trust-region-reflective
	%All-forwards fit obtained: 	-5.0653      -3.9939	-3.7181		-2.7503		-3.0933		-3.9808		1.669													 with norm 0.35161 in 110.995 seconds using levenberg-marquardt*
	%All fit obtained: 			-5.0711      -3.583		-3.2818		-3.6145		-4.9342		-3.138		-3.7425     -3.2157     -2.9387     -5.9996 			 with norm 0.86563 in 485.6 seconds using trust-region-reflective
	%All fit obtained: 			-5.0216      -3.5615	-3.2697		-3.6269		-4.936		-3.1299		-3.922      -3.2063     -2.9664     -7.014 				 with norm 0.86588 in 459.6819 seconds using levenberg-marquardt
	%All-diffusion fit obtained: -5.1	     -3.5917	-3.2835		-3.6365		-4.9349		-3.1374		-3.666      -3.2155     -2.9604     -5.9995    -6.9209   with norm 0.86557 in 705.7446 seconds using trust-region-reflective
	%All-diffusion fit obtained: -5.0239      -3.5641	-3.271		-3.7369		-4.9272		-3.1277		-3.7986     -3.2121     -3.1069      -7.12     -3.9932   with norm 0.86585 in 534.7464 seconds using levenberg-marquardt
	%
	%Initial conditions: -6 -6 -6... -2 -2 -2...
	%
	%Parameters:
	%---------------------
	%Base fit obtained:			-3.3139		-3.0083		-3.6464		-3.3433																						 with norm 0.45566 in 204.3315 seconds using trust-region-reflective
	%Base fit obtained:			-3.1385		-2.9245		-3.6407		-5.7897																						 with norm 0.45618 in 120.15 seconds using levenberg-marquardt*
	%Diffusion fit obtained:		-3.1344		-2.973		-3.6561		-3.4464		-7.6648																			 with norm 0.4648 in 245.7179 seconds using trust-region-reflective
	%Diffusion fit obtained:		-2.7091		-2.2902		-3.6965		-4.4702		-8.1513																			 with norm 0.46918 in 225.0802 seconds using levenberg-marquardt
	%All-forwards fit obtained: 	-1.7133		-3.0183		-1.7799		-1.0357		-0.053597	-0.92305	-1.4226													 with norm -1.0578 in 229.3659 seconds using trust-region-reflective
	%All-forwards fit obtained: 	-2.5833		-3.51455	-3.32525	-2.32137	-3.11419	-3.42393	-13.4732												 with norm 0.67836 in 284.6964 seconds using levenberg-marquardt*
	%All fit obtained:			-2			-2.0001		-2		  	-2	 		-2.0001		-6			-6			-6		  	-6		  	-6					 with norm -0.47303 in 120.488 seconds using trust-region-reflective
	%All fit obtained: 			-2.1444		-1.8747		-1.9949	  	-2.0009		-2.279		-5.9789		-6			-5.9931	  	-5.9997		-6					 with norm -0.28761 in 201.7314 seconds using levenberg-marquardt*
	%All-diffusion fit obtained: -2.0759		-1.9376		-2.0002		-1.9997		-2.0051		-1.9127		-5.1365		-4.1552		-6			-5.4805		-5.7477  with norm -0.22214 in 385.1773 seconds using trust-region-reflective
	%All-diffusion fit obtained: -1.40996	33.3311		-37.0317	-2.01164	 2.32856	-2.89004	-6.00026	2.99014		6.67349	  	4.61599		-10.1654 with norm 0.21565 in 317.0109 seconds using levenberg-marquardt
	%
	% * = some solutions had errors
	
	experiment = direct_exp4o_lsq;
	experiment_diff = direct_exp4o_diff_lsq;
	filename_tbim = '';
	filename_tbid = '';
	plotfit = 0;
	verbose = 0;
	method1 = 'trust-region-reflective';
	method2 = 'levenberg-marquardt';
	
	%Compare base model. 4 free parameters
	x0 = [-4 -4 -4 -4];
	tstart = tic;
	[fits resnorm] = fit_exp4(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'Base');
	tstart = tic;
	[fits resnorm] = fit_exp4(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
	print_results(fits, resnorm, tstart, method2, 'Base');
	
	%Compare base model with diffusion coefficient. 5 free parameters
	x0 = [-4 -4 -4 -4 -4];
	tstart = tic;
	[fits resnorm] = fit_exp4_diff(experiment_diff, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'Diffusion');
	tstart = tic;
	[fits resnorm] = fit_exp4_diff(experiment_diff, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
	print_results(fits, resnorm, tstart, method2, 'Diffusion');
	
	%Compare 'all' forward model. 7 free parameters
	x0 = [-4 -4 -4 -4 -4 -4 -4];
	tstart = tic;
	[fits resnorm] = fit_exp4_allforward(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'All-forwards');
	tstart = tic;
	[fits resnorm] = fit_exp4_allforward(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
	print_results(fits, resnorm, tstart, method2, 'All-forwards');
	
	%Compare 'all' model. 10 free parameters
	x0 = [-4 -4 -4 -4 -4 -4 -4 -4 -4 -4];
	tstart = tic;
	[fits resnorm] = fit_exp4_all(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'All');
	tstart = tic;
	[fits resnorm] = fit_exp4_all(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
 	print_results(fits, resnorm, tstart, method2, 'All');
	
	%Compare 'all' model with diffusion coefficient. 11 free parameters
	x0 = [-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4];
	tstart = tic;
	[fits resnorm] = fit_exp4_all_diff(experiment_diff, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'All-diffusion');
	tstart = tic;
	[fits resnorm] = fit_exp4_all_diff(experiment_diff, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
	print_results(fits, resnorm, tstart, method2, 'All-diffusion');

	%Run same tests but with less realistic starting values...
	%Compare base model. 4 free parameters
	x0 = [-2 -2 -6 -6];
	tstart = tic;
	[fits resnorm] = fit_exp4(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'Base');
	tstart = tic;
	[fits resnorm] = fit_exp4(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
	print_results(fits, resnorm, tstart, method2, 'Base');
	
	%Compare base model with diffusion coefficient. 5 free parameters
	x0 = [-2 -2 -2 -6 -6];
	tstart = tic;
	[fits resnorm] = fit_exp4_diff(experiment_diff, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'Diffusion');
	tstart = tic;
	[fits resnorm] = fit_exp4_diff(experiment_diff, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
	print_results(fits, resnorm, tstart, method2, 'Diffusion');
	
	%Compare 'all' forward model. 7 free parameters
	x0 = [-2 -2 -2 -2 -6 -6 -6];
	tstart = tic;
	[fits resnorm] = fit_exp4_allforward(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'All-forwards');
	tstart = tic;
	[fits resnorm] = fit_exp4_allforward(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
	print_results(fits, resnorm, tstart, method2, 'All-forwards');
	
	%Compare 'all' model. 10 free parameters
	x0 = [-2 -2 -2 -2 -2 -6 -6 -6 -6 -6];
	tstart = tic;
	[fits resnorm] = fit_exp4_all(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'All');
	tstart = tic;
	[fits resnorm] = fit_exp4_all(experiment, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
 	print_results(fits, resnorm, tstart, method2, 'All');
	
	%Compare 'all' model with diffusion coefficient. 11 free parameters
	x0 = [-2 -2 -2 -2 -2 -2 -6 -6 -6 -6 -6];
	tstart = tic;
	[fits resnorm] = fit_exp4_all_diff(experiment_diff, x0, filename_tbim, filename_tbid, plotfit, verbose, method1);
	print_results(fits, resnorm, tstart, method1, 'All-diffusion');
	tstart = tic;
	[fits resnorm] = fit_exp4_all_diff(experiment_diff, x0, filename_tbim, filename_tbid, plotfit, verbose, method2);
	print_results(fits, resnorm, tstart, method2, 'All-diffusion');
end

function print_results(fits, resnorm, tstart, method, model)
	disp([model ' fit obtained: ' num2str(fits) ' with norm ' num2str(resnorm) ' in ' num2str(toc(tstart)) ' seconds using ' method]);
end