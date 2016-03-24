function fit_mcmc(data_out, nruns)
	% fit_mcmc	Run MCMC using http://helios.fmi.fi/~lainema/mcmc/ to generate
	%			estimate of posterior
	%
	% Usage:
	%			fit_mcmc(nruns, data_out)
	%
	% Input:
	%			nruns
	%			data_out = name of output file to write data to

	%Setup params
	%A0 = 0.02090; B0 = A0/3; C0 = 0; D0 = 0; E0 = 0;
	%data.y0 = [A0;B0;C0;D0;E0];

	%k00 = [15,1.5,0.3]';
	%[k0,ss0] = fminsearch(@himmelss,k00,[],data)
	%mse = ss0/(length(data.ydata)-4);

	%   {  {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
	%      {'par2',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
	%      ... }

	if (nargin < 2) nruns = 5000; end 
	assert(nargin > 0, 'Please specify output file name')

	nW = 300;
	nPars = 25;
	minit = zeros(nPars, nW);

	%test run
	%nruns = 10;
	%nW = 10;

	[obs, obs_std] = observations_exp46();

	ssfun   = @(k) likelihood(k, obs, obs_std);
	priorfun = @(k) prior(k);
	logPfuns = {priorfun ssfun};	
	
	%make a set of starting points for the entire ensemble of walkers 
	for idx = 1:nW
		minit(:,idx) = sampleprior; 
	end
	
	%Apply the MCMC hammer 
	[models,logP]=gwmcmc(minit,logPfuns,nW*nruns); 
	%models(:,:,1:floor(size(models,3)*.2))=[]; %remove 20% as burn-in 
	%models=models(:,:)'; %reshape matrix to collapse the ensemble member dimension 

	%for j = 1:nPars
	%	for i = j:nPars
	%		subplot(nPars*nPars, j, i)
	%		scatter(models(:,i),models(:,j))
	%	end
	%end

	%scatter(models(:,1),models(:,2)) 
	%prctile(models,[5 50 95]) 

	save(data_out, 'models', 'logP');

	%f_tbid = run_simulation(model, 'tBid', T_p_exp14('tBid'), 0:(181*60), fits(tbid_params));
	%f_tbim = run_simulation(model, 'tBim', T_p_exp14('tBim'), 0:(181*60), fits(tbim_params));
	%[r2 bakmcl1r2bid tbidmcl1r2 bakmcl1r2bim tbimmcl1r2] = compute_r2_exp14(observations, resnorm, f_tbid, f_tbim, time_points);
	%if (nargin < 4) | (plotfit == 1)
	%	plot_fits('tBim', model, fits(tbim_params), observations(1:7,:), filename_tbim);
	%	plot_fits('tBid', model, fits(tbid_params), observations(8:14,:), filename_tbid);
	%end
end

function plot_fits(activator, model, x, observations, filename)
	t = 0:(180*60);
	time_points = [0 10 20 30 40 50 60 70 80 90 120 150];
 	[f, h1, h2, h3, h4, h5] = run_simulation(model, activator, T_p_exp14(activator), t, x);
 	
	colororder = get(0, 'DefaultAxesColorOrder');
	figure(h1);
	xlim([0 180]);
	%Bak:Mcl-1
	line(time_points, observations(3,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(3,:));
	%Bak multimer
	%line(time_points, observations(6,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(5,:));
	saveplot(h1, [filename '-bak.eps']);

	figure(h2);
	xlim([0 180]);
	line(time_points, observations(7,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(7,:));
	saveplot(h2, [filename '-activator.eps']);

	figure(h3);
	xlim([0 180]);
	%Bak:Mcl-1
	line(time_points, observations(3,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(3,:));
	%tBid:Mcl-1
	line(time_points, observations(7,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(7,:));
	saveplot(h3, [filename '-mcl-1.eps']);

	figure(h4);
	xlim([0 180]);
	line(time_points, observations(1,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(9,:));
	line(time_points, observations(2,:), 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 3, 'Color', colororder(10,:));
	saveplot(h4, [filename '-cytc.eps']);
	
	figure(h5);
	xlim([0 180]);
	saveplot(h5, [filename '-bclxl.eps']);
end
