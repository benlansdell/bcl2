function R = robust_analysis_bentele_rndparam_6(n, basename)
	% robust_analysis_bentele_rndparam_6		Function to analyse the sensitivity of model defined by specified
	%						kinetic parameters. Picks parameters randomly and collects statistics
	%						on sensitivity.  Uses 'global' method from Bentele 04. 
	%						Also compute the Hessian matrix for the goodness of fit, 
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
	%						R.hessian = matrix of the 'mean Hessian'
	% 
	% Examples:
	%						R = robust_analysis_bentele_rndparam_6_heatmap(2, './images/robust_analysis_bentele_rndparam_6');

	%Change parameter by 10%
	dk = 0.99;
	activator = 'tBid';
	model = 'random_parameters_both';

	variance = 500;
	observations = observations_exp14_multi();
	obs_tbim = observations(1:7,:);
	obs_tbid = observations(8:14,:);
	if strcmp(activator, 'tBid')
		obs = obs_tbid;
	else
		obs = obs_tbim;
	end

	%Direct activation model
	x = [1 1 2 2 3 3 4 4 5 6  6 7];
	y = [1 2 1 2 1 2 1 2 1 1  2 1];
	z = [1 2 3 4 5 6 7 8 9 10 11 12 13];
	[k, ic, names] = load_model6(activator, model);
	t = (-50*60):(120*60);
	f = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
	S = zeros(n,13,12);
	Sweight = zeros(n,13,12);
	SSR = zeros(n,1);

	matlabpool 4;
	parfor m = 1:n
		[k, ic, names] = load_model6(activator, model);
		f = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
		SSR(m) = ss_resid(k, obs, activator, model);
		display(['Run ' num2str(m) '. SSR=' num2str(SSR(m))]);
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
		%Weigh all entries by exp(-SSR)
		Sweight(m,:,:) = S(m,:,:)*exp(-SSR(m)/variance);
	end
	%Normalize
	%Sweight = Sweight/exp(-sum(SSR));
	matlabpool close;

	R = struct('mu', reshape(mean(S,1),13,12), 'sigma', reshape(std(S,0,1),13,12));
	R.muweight = reshape(mean(Sweight,1),13,12);
	R.sigmaweight = reshape(std(Sweight,0,1),13,12);
	R.SSR = SSR;
	colororder = get(0, 'DefaultAxesColorOrder');

	paramnames = {'k1a', 'k1d', 'k2a', 'k2d', 'k3a', 'k3d', 'k4a', 'k4d', 'k5a', 'k6a', 'k6d', 'k7a'};
	imagesc(log(R.mu+0.001));
	set(gca,'XTick',1:length(paramnames));
	set(gca,'XTickLabel',paramnames);
	set(gca,'YTick',1:length(names));
	set(gca,'YTickLabel',names);
	ylabel('Species');
	xlabel('Parameters');
	colorbar
	saveplot(gcf, [basename '_mu.eps']);

	imagesc(log(R.sigma+0.001));
	set(gca,'XTick',1:length(paramnames));
	set(gca,'XTickLabel',paramnames);
	set(gca,'YTick',1:length(names));
	set(gca,'YTickLabel',names);
	ylabel('Species');
	xlabel('Parameters');
	colorbar
	saveplot(gcf, [basename '_sigma.eps']);
      
	imagesc(log(R.sigmaweight+0.001));
	set(gca,'XTick',1:length(paramnames));
	set(gca,'XTickLabel',paramnames);
	set(gca,'YTick',1:length(names));
	set(gca,'YTickLabel',names);
	ylabel('Species');
	xlabel('Parameters');
	colorbar
	saveplot(gcf, [basename '_sigmaweight.eps']);

	imagesc(log(R.muweight+0.001));
	set(gca,'XTick',1:length(paramnames));
	set(gca,'XTickLabel',paramnames);
	set(gca,'YTick',1:length(names));
	set(gca,'YTickLabel',names);
	ylabel('Species');
	xlabel('Parameters');
	colorbar
	saveplot(gcf, [basename '_muweight.eps']);

	save('./scripts/robust_analysis_bentele_rndparam_6_Smatrix.mat', 'S', 'Sweight');
end

function ss = ss_resid(k, observations, activator, model)
	%x is parameters vector, observations is observations matrix, model is not too important, and 
	resid = [solve_ode(activator, model, k)] - observations;
	ss = sum(sum(resid.^2));
end

function f_scaled = solve_ode(activator, model, k)	
	t = 0:(181*60);
	time_points = 60*[0 10 20 30 40 50 60 70 80 90 120 150] + 1;
	f = run_simulation(model, activator, T_p_exp14(activator), t, k);
	
	cytcc = f(9, time_points);
	cytcm = f(5, time_points);	
	bakmcl1 = f(8, time_points);
	tbimmcl1 = f(7, time_points);
	bak = f(1, time_points);
	abak = f(6, time_points);
	dimbak = f(4, time_points);
	mulbak = f(10, time_points);
	
	f_scaled = [cytcc; cytcm; bakmcl1; abak; dimbak; mulbak; tbimmcl1];
	
end

function f = run_simulation(model, activator, T_p, t, k)
	[x, ic, names] = load_model14(activator, model);
	f = bcl2model(k, ic, names, t, T_p, activator);
end