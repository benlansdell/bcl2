function diffusion_exp6o(activator)
	% diffusion_exp6o	Plot solutions to Exp 6 model for a range of diffusive constants. Uses dimensional variables.
	%			Specify which set of parameters to use.
	%
	% Usage:
	%			diffusion_exp6o(param_mode, activator)
	%
	% Input:
	%			activator(optional) = one of 'tBid' or 'tBim'. Default 'tBid'
	%
	% Example(s):
	%			diffusion_exp6o('tBid');

	%tBIM
	obs_tbim = [0.043105305	0.000693856	0.09873329	0.01417368	0.282553204	0.999623494	1.622121962	2.516921961	2.714091981	3.270253467	3.414746138	3.554614254	3.586706398	4.143876929	4.508449773; % cyt c super
			4.956894695	4.999306144	4.90126671	4.98582632	4.717446796	4.000376506	3.377878038	2.483078039	2.285908019	1.729746533	1.585253862	1.445385746	1.413293602	0.856123071	0.491550227; % cyt c mitos
			0.781	0.846	2.577	10.143	4.027	2.564	1.600	2.899	0.073	0.775	0.748	0.642	0.991	0.293	0.358; %Mcl-1:Bak
			6.859273759	11.04637059	13.19871582	20	20	20	10.71630187	19.55938416	20	20	20	20	20	20	19.64216667]; %tBIM:Mcl-1

	%tBid
	obs_tbid = [0.003840869	0.010521271	0.006294496	0.018543081	0.028324882	0.071314757	0.097875765	0.233443547	0.314826241	0.376005415	0.450318083	0.455888639	0.459954039	0.899681844	1.057177511; %cyt c super
			4.996159131	4.989478729	4.993705504	4.981456919	4.971675118	4.928685243	4.902124235	4.766556453	4.685173759	4.623994585	4.549681917	4.544111361	4.540045961	4.100318156	3.942822489; %cyt c mito
			0.649	0.764	0.823	3.424	5.260	4.368	4.191	6.917	6.452	6.439	6.871	5.558	2.588	5.847	6.769; %Mcl-1:Bak
			0.99561264	1.076035389	2.042208005	2.70293412	3.890831606	3.628270574	5.628992636	11.52109612	9.902594541	8.56284851	12.55884448	13.87699783	9.534669212	15.40522999	13.23116667]; %tBid:Mcl-1

	model = 'reversible_exp6_diff';
	if (nargin < 1)
		activator = 'tBid';
	end
	D = [100 10^(-8) 10^(-8.5) 10^(-9) 10^(-9.5) 10^(-10)];
	[k, ic, names] = load_model6(activator, model);
	T_p = T_p_exp6(activator);
	t = 0:10800;

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
	time_points = [0 10 20 30 40 50 60 70 80 90 100 110 120 150 180];
	prev_ColorOrder = get(0, 'DefaultAxesColorOrder');
	my_ColorOrder = [1.00000   0.00000   0.00000;
				     0.80000   0.00000   0.20000;
				     0.60000   0.00000   0.40000;
	  			     0.40000   0.00000   0.60000;
	  			     0.20000   0.00000   0.80000;
	  				 0.00000   0.00000   1.00000;];
	%set(0,'DefaultAxesColorOrder',my_ColorOrder);

	plot(t/60, bakmcl1)
	if strcmp(activator, 'tBid')
		line(time_points, obs_tbid(3,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', 'b');
	else
		line(time_points, obs_tbim(3,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', 'b');	
	end
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 15]);
	%legend('D = \infty', 'D = 10^{-8}', 'D = 10^{-8.5}', 'D = 10^{-9}', 'D = 10^{-9.5}', 'D = 10^{-10}', 'Location', 'NorthEast');

	filename = ['./images/demos/diffusion_exp6o_bakmcl1-' activator '.eps'];
	saveplot(gcf, filename);

	plot(t/60, tbidmcl1)
	if strcmp(activator, 'tBid')
		line(time_points, obs_tbid(4,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', 'b');
	else
		line(time_points, obs_tbim(4,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', 'b');	
	end
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 23]);
	%legend('D = \infty', 'D = 10^{-8}', 'D = 10^{-8.5}', 'D = 10^{-9}', 'D = 10^{-9.5}', 'D = 10^{-10}', 'Location', 'NorthWest');

	filename = ['./images/demos/diffusion_exp6o_tbidmcl1-' activator '.eps'];
	saveplot(gcf, filename);

	plot(t/60, dimbak)
	xlabel('time (min)');
	ylabel('concentration (nM)');
	ylim([0 3.5]);
	%legend('D = \infty', 'D = 10^{-8}', 'D = 10^{-8.5}', 'D = 10^{-9}', 'D = 10^{-9.5}', 'D = 10^{-10}', 'Location', 'NorthWest');

	filename = ['./images/demos/diffusion_exp6o_bak-' activator '.eps'];
	saveplot(gcf, filename);

	set(0, 'DefaultAxesColorOrder', prev_ColorOrder);
