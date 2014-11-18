function diffusion_exp3o(param_mode, T0)
	% diffusion_exp3o	Plot solutions to Exp 3o model for a range of diffusive constants. Uses dimensional variables.
	%			Specify which set of parameters to use. Only plot for tBid case
	%
	% Usage:
	%			diffusion_exp3o(param_mode)
	%
	% Input:
	%			param_mode = one of 'all' -- 
	%					 or 'naive' --
	%					 or 'ratio' --
	%			T0 = initial concentration of tBid
	%
	% Example(s):
	%			diffusion_exp3o('all', 10);

	experiment = direct_exp3o_diff_lsq;
	activator = 'tBid';
	D = [100 10^(-8) 10^(-8.5) 10^(-9) 10^(-9.5) 10^(-10)];

	[k, ic, t] = experiment.setup_params(activator);
	ic(2) = T0;
	if strcmp(param_mode, 'all')
		k(1,1) = 10^(-3.583);
		k(2,1) = 10^(-3.2818);
		k(2,2) = 10^(-3.2157);
		k(3,1) = 10^(-4.9342);
		k(3,2) = 10^(-5.9996);
		k(4,1) = 10^(-3.138);
		k(5,1) = 10^(-3.7425);
	elseif strcmp(param_mode, 'naive')
		k(1,1) = 10^(-3.0500);
		k(2,1) = 10^(-2.8633);
		k(2,2) = 10^(-3.0980);
		k(3,1) = 10^(-2.0915);
		k(3,2) = 10^(-1.1238);
		k(4,1) = 10^(-3.8100);
		k(5,1) = 10^(-5.9200);
	elseif strcmp(param_mode, 'ratio')
		k(1,1) = 10^(-3.5100);
		k(2,1) = 10^(-3.4600);
		k(2,2) = 10^(-3.2100);
		k(3,1) = 10^(-4.8400);
		k(3,2) = 10^(-3.6000);
		k(4,1) = 10^(-3.5700);
		k(5,1) = 10^(-4.1100);
	else
		throw(MException('ParamError:invalidParameters', 'Invalid param_mode: must be one of all, naive or ratio'));
	end

	bakmcl1 = zeros(length(D), length(t));
	bak = zeros(length(D), length(t));
	tbidmcl1 = zeros(length(D), length(t));
	
	%Solve for different values of D
	for idx = 1:length(D)
		f1 = experiment.direct_exp3(k, ic, D(idx), t, activator, 0);
		bakmcl1(idx,:) = f1(8,:);
		bak(idx,:) = f1(1,:);
		tbidmcl1(idx,:) = f1(7,:);
	end
	
	figure
	prev_ColorOrder = get(0, 'DefaultAxesColorOrder');
	my_ColorOrder = [0.00000   0.00000   0.00000;
				     0.00000   0.00000   0.20000;
				     0.00000   0.00000   0.40000;
	  			     0.00000   0.00000   0.60000;
	  			     0.00000   0.00000   0.80000;
	  				 0.00000   0.00000   1.00000;];
	set(0,'DefaultAxesColorOrder',my_ColorOrder);

	plot(t/60, bakmcl1)
	xlabel('time (min)');
	ylabel('Bak:Mcl-1 concentration (nM)');
	ylim([0 15]);
	legend('D = \infty', 'D = 10^{-8}', 'D = 10^{-8.5}', 'D = 10^{-9}', 'D = 10^{-9.5}', 'D = 10^{-10}', 'Location', 'NorthEast');

	filename = ['./images/demos/diffusion_exp3o_bakmcl1-' param_mode '.eps'];
	set(gcf, 'paperunits', 'inches');
	set(gcf, 'papersize', [6 4]);
	set(gcf, 'paperposition', [0 0 6 4]);
	print(gcf,  filename, '-depsc');

	plot(t/60, tbidmcl1)
	xlabel('time (min)');
	ylabel('tBid:Mcl-1 concentration (nM)');
	ylim([0 23]);
	legend('D = \infty', 'D = 10^{-8}', 'D = 10^{-8.5}', 'D = 10^{-9}', 'D = 10^{-9.5}', 'D = 10^{-10}', 'Location', 'NorthWest');

	filename = ['./images/demos/diffusion_exp3o_tbidmcl1-' param_mode '.eps'];
	set(gcf, 'paperunits', 'inches');
	set(gcf, 'papersize', [6 4]);
	set(gcf, 'paperposition', [0 0 6 4]);
	print(gcf,  filename, '-depsc');

	plot(t/60, bak)
	xlabel('time (min)');
	ylabel('Bak concentration (nM)');
	ylim([0 15]);
	legend('D = \infty', 'D = 10^{-8}', 'D = 10^{-8.5}', 'D = 10^{-9}', 'D = 10^{-9.5}', 'D = 10^{-10}', 'Location', 'NorthEast');

	filename = ['./images/demos/diffusion_exp3o_bak-' param_mode '.eps'];
	set(gcf, 'paperunits', 'inches');
	set(gcf, 'papersize', [6 4]);
	set(gcf, 'paperposition', [0 0 6 4]);
	print(gcf,  filename, '-depsc');

	set(0, 'DefaultAxesColorOrder', prev_ColorOrder);
