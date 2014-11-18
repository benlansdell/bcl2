function R = robust_analysis_bentele_14_heatmap(model, activator, filename)
	% robust_analysis_bentele		Function to analyse the sensitivity of model defined by specified
	%						kinetic parameters. Use 'global' method from Bentele 04. 
	%						Compute also Hessian matrix which can be used to measure
	%						sloppiness of model near given parameters
	%
	% Usage:
	%
	%						R = robust_analysis_bentele_14_heatmap(model, activator, filename)
	%
	% Input:
	%						k = 9x2 matrix of kinetic parameters
	%						ic = 13x1 vector of initial concentrations
	%						filename = output EPS file for 2d bar plot
	%
	% Output:
	%						R = structure with the following components:
	%						R.S = Sensitivity matrix
	%						R.H = Hessian matrix
	% 
	% Examples:
	%						R = robust_analysis_bentele_14_heatmap('reversible_exp14_both', 'tBid', './images/robust_analysis_bentele_14heatmap_tbid.eps');

	dk = 1.01;
	logdk = log10(dk);
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
	[k, ic, names] = load_model14(activator, model);
	t = (-50*60):(120*60);
	f = bcl2model(k, ic, names, t, T_p_exp14(activator), activator);
	SSR = ss_resid(k, obs, activator, model);
	S = zeros(13,12);
	chi2 = zeros(12,12);
	%For each parameter
	for i=1:12
		i
		display('Monitoring change in species sensitivity')
		%For each species compute sensitivity
		for j=1:13
			j
			k_ = k;
			k_(x(i), y(i)) = k(x(i), y(i))*dk;
			f_ = bcl2model(k_, ic, names, t, T_p_exp14(activator), activator);
			c = f(z(j),:);
			c_ = f_(z(j),:);
			C = sum(c, 2);
			C_ = sum(abs(c_-c),2);
			S(j,i) = -log10(k(x(i),y(i)))*C_/C/abs(logdk);
		end
		%For each other parameter
		display('Computing Hessian entries')
		%Symmetric matrix
		for j = i:12			
			j
			ki = k;
			kj = k;
			kij = k;
			ki(x(i), y(i)) = k(x(i), y(i))*dk;
			kj(x(j), y(j)) = k(x(j), y(j))*dk;
			kij(x(i), y(i)) = k(x(i), y(i))*dk;
			kij(x(j), y(j)) = k(x(j), y(j))*dk;
			%Compute SSR for perturbed models
			SSRi = ss_resid(ki, obs, activator, model);
			SSRj = ss_resid(kj, obs, activator, model);
			SSRij = ss_resid(kij, obs, activator, model);
			%Compute Hessian entry
			%chi2(i,j) = ki(x(i),y(i))*kj(x(j),y(j))*(SSRij-SSRi-SSRj+SSR)/logdk^2;
			chi2(i,j) = (SSRij-SSRi-SSRj+SSR)/logdk^2;
			chi2(j,i) = chi2(i,j);
		end
	end
	R = struct('S', S);
	R.H = chi2;
	
	colororder = get(0, 'DefaultAxesColorOrder');
	paramnames = {'k1a', 'k1d', 'k2a', 'k2d', 'k3a', 'k3d', 'k4a', 'k4d', 'k5a', 'k6a', 'k6d', 'k7a'};
	imagesc(log(S+0.001));
	set(gca,'XTick',1:size(S,2));
	set(gca,'XTickLabel',paramnames);
	set(gca,'YTick',1:size(S,1));
	set(gca,'YTickLabel',names);
	xlabel('Parameters');
	ylabel('Species');
	colorbar
	saveplot(gcf, filename);
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