function diffusion_exp4o(activator)
	% diffusion_exp4o	Plot solutions to Exp 4 model for a range of diffusive constants. Uses dimensional variables.
	%			Specify which set of parameters to use.
	%
	% Usage:
	%			diffusion_exp4o(param_mode, activator)
	%
	% Input:
	%			activator(optional) = one of 'tBid' or 'tBim'. Default 'tBid'
	%
	% Example(s):
	%			diffusion_exp4o('tBid');

	obs = [ 0.132578077	0.152173913	0.16442131	0.184935701	0.253521127	0.330679731	0.282914881	0.345988977	0.407225964;	%tBid: Cyt C S/N release
			3.123086344	3.919167177	4.378444581	3.85793019	3.919167177	4.500918555	4.654011023	4.654011023	3.704837722;	%tBid: Cyt C mitos
			0.917647059	0.952941176	1.561764706	5.426470588	9			8.241176471	9			7.173529412	5.902941176;	%tBid: Bak*:Mcl-1
			1			1			1			1			1			1			1			1			1;				%tBid: Bak* free
			1			1			1			1			1			1			1			1			1;				%tBid: Bak* dimer
			0.162781955	0.186842105	0.201879699	0.227067669	0.311278195	0.406015038	0.347368421	0.42481203	0.5;			%tBid: Bak* multimer -- same form as Cyt C S/N release but rescaled
			0.338928571	0.314166667	0.598154762	0.858928571	1.516666667	1.973214286	6.360714286	15.70833333	13];			%tBid: tBid:Mcl-1

	model = 'reversible_exp4_diff';
	if (nargin < 1)
		activator = 'tBid';
	end
	D = [100 10^(-8) 10^(-8.5) 10^(-9) 10^(-9.5) 10^(-10)];
	[k, ic, names] = load_model(activator, model);
	T_p = T_p_exp4(activator);
	t = 0:5400;

	bakmcl1 = zeros(length(D), length(t));
	bak = zeros(length(D), length(t));
	tbidmcl1 = zeros(length(D), length(t));
	dimbak = zeros(length(D), length(t));
	
	%Solve for different values of D
	for idx = 1:length(D)
		conc = bcl2model_diff(k, ic, D(idx), names, t, T_p, activator);
		dimbak(idx,:) = conc(4,:);
		bakmcl1(idx,:) = conc(8,:);
		bak(idx,:) = conc(1,:);
		tbidmcl1(idx,:) = conc(7,:);
	end

	figure
	time_points = [0 10 20 30 40 50 60 70 80];
	prev_ColorOrder = get(0, 'DefaultAxesColorOrder');
	my_ColorOrder = [0.00000   0.00000   0.00000;
				     0.00000   0.00000   0.20000;
				     0.00000   0.00000   0.40000;
	  			     0.00000   0.00000   0.60000;
	  			     0.00000   0.00000   0.80000;
	  				 0.00000   0.00000   1.00000;];
	set(0,'DefaultAxesColorOrder',my_ColorOrder);

	plot(t/60, bakmcl1)
	if strcmp(activator, 'tBid')
		line(time_points, obs(3,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', 'b');
	end
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 15]);
	legend('D = \infty', 'D = 10^{-8}', 'D = 10^{-8.5}', 'D = 10^{-9}', 'D = 10^{-9.5}', 'D = 10^{-10}', 'Location', 'NorthEast');

	filename = ['./images/demos/diffusion_exp4o_bakmcl1-' activator '.eps'];
	saveplot(gcf, filename);

	plot(t/60, tbidmcl1)
	if strcmp(activator, 'tBid')
		line(time_points, obs(7,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', 'b');
	end
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 23]);
	legend('D = \infty', 'D = 10^{-8}', 'D = 10^{-8.5}', 'D = 10^{-9}', 'D = 10^{-9.5}', 'D = 10^{-10}', 'Location', 'NorthWest');

	filename = ['./images/demos/diffusion_exp4o_tbidmcl1-' activator '.eps'];
	saveplot(gcf, filename);

	plot(t/60, dimbak)
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 3.5]);
	legend('D = \infty', 'D = 10^{-8}', 'D = 10^{-8.5}', 'D = 10^{-9}', 'D = 10^{-9.5}', 'D = 10^{-10}', 'Location', 'NorthWest');

	filename = ['./images/demos/diffusion_exp4o_bak-' activator '.eps'];
	saveplot(gcf, filename);

	set(0, 'DefaultAxesColorOrder', prev_ColorOrder);