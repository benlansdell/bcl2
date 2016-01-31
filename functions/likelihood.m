function logl = likelihood(x, observations, observations_std)
	%ssfun(par,data)

	% likelihood	Compute log-likelihood for a given set of parameters by simulating
	%				mass-action kinetic model of protein-protein interactions. x is typically
	%				sampled from MCMC routine
	%
	% Usage:
	%			logl = likelihood(x)
	%
	% Input:
	%			x = parameters for k1a_tbim, k1a_tbid, k2a, k3a_tbim, k4a, k5a, k6a, k8a, k9a_tbim, k2d, k3d_tbim, k4d, k6d, k8d, k9d_tbim
	%
	% Output:
	%			logl = log-likelihood

%	   1 		      3         4              6    7    8    9   10        11             13   14
%	    		 2    3                   5    6    7    8    9   10                  12   13   14
%k1a_tbim, k1a_tbid, k2a, k3a_tbim, k3a_tbid, k4a, k5a, k6a, k7a, k8a, k9a_tbim, k9a_tbid, k1d, k2d,
%      15             17   18   19   20   21        22        
%                16   17   18   19   20   21                  23
%k3d_tbim, k3d_tbid, k4d, k5d, k6d, k7d, k8d, k9d_tbim, k9d_tbid

	model = 'reversible_exp46';
	tbim_params = [1 3 4 6 7 8 9 10 11 13 14 15 17 18 19 20 21 22 24 25];
	tbid_params = [2 3 5 6 7 8 9 10 12 13 14 16 17 18 19 20 21 23 24 25];

	[obs, obs_std] = scale(observations, observations_std, x, model);
	ssr = ss_resid(x, obs, obs_std, model, tbim_params, tbid_params);
	logl = -ssr/2;
end

function [obs, obs_std] = scale(observations, observations_std, x, model)
	%ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];
	[k, ic, names] = load_model46('tBim', model);
	B_0 = x(24)
	M_0 = ic(3);
	C_0 = ic(5);
	%Scale the Bak limited observations
	baklimited = [2, 3, 6, 7];
	cytclimited = [1, 5];
	mcllimited = [4, 8];

	obs = observations;
	obs_std = observations_std;

	obs(baklimited, :) = observations(baklimited,:)/100*B_0;
	obs(cytclimited, :) = observations(cytclimited,:)/100*C_0;
	obs(mcllimited, :) = observations(mcllimited,:)/100*M_0;
	obs_std(baklimited, :) = observations_std(baklimited,:)/100*B_0;
	obs_std(cytclimited, :) = observations_std(cytclimited,:)/100*C_0;
	obs_std(mcllimited, :) = observations_std(mcllimited,:)/100*M_0;

	%	6.3, 2.6, 4.5, 6.5, 15.3, 20.5, 26.9, 37.4, 64.2, 96.0; 							%BidBim: CytC 		(CytC limited. x 5/100)
	%	9.9, 8.4, 11.8, 36.9, 95.6, 97.1, 95.2, 72.6, 65.9, 26.1; 							%Bak:Mcl1 			(Bak limited. x B_0/100)
	%	100.0, 92.0, 83.3, 48.2, 5.0, 2.3, 2.4, 2.1, 2.7, 2.7; 								%Inactive Bak 		(Bak limited. x B_0/100)
	%	0.03*100/35 0.27*100/35 27*100/35	100	  100	100    100 	100	  100	100; 		%tBim:Mcl-1 		(Mcl-1 limited. x 35/100)
	%	1.0, 4.4, 5.6, 11.7, 14.1, 22.9, 24.0, 23.3, 21.9, 31.6; 							%Bid: CytC release	(CytC limited. x 5/100)
	%	11.8, 9.9, 19.5, 78.2, 91.8, 78.6, 86.1, 92.1, 67.4, 63.9; 							%Bak:Mcl1			(Bak limited. x B_0/100)
	%	100.0, 94.9, 64.1, 11.1, 2.9, 1.5, 1.2, 1.1, 1.2, 1.0;								%Inactive Bak		(Bak limited. x B_0/100)
	%	0.03*100/35 0.27*100/35 27*100/35	100	  100	100    100 	100	  100	100]; 		%tBid:Mcl-1			(Mcl-1 limited. x 35/100)
end

function ssr = ss_resid(x, observations, obs_std, model, tbim_params, tbid_params)
	res = residuals(x, observations, model, tbim_params, tbid_params);
	%res
	obs_std 
	%res.^2./(obs_std.^2)
	ssr = sum(sum(res.^2./(obs_std.^2)));
end

function res = residuals(x, observations, model, tbim_params, tbid_params)
	%Simulate both tBim and tBid experiments, return log-likelihood from both
	%TODO (eventually): add extra experiments with out Mcl1
	res = [solve_ode('tBim', model, x(tbim_params));
	       solve_ode('tBid', model, x(tbid_params))] - observations;
end

function f_scaled = solve_ode(activator, model, x)	
	t = (-50*60):(181*60);
	time_points = (50*60) + 60*[0 10 20 30 40 50 60 90 120 150] + 1;
	f = run_simulation(model, activator, T_p_exp46(activator), t, x);	
	cytcc = f(9, time_points);
	bakmcl1 = f(8, time_points);
	tbimmcl1 = f(7, time_points);
	bak = f(1, time_points);
	abak = f(6, time_points);
	dimbak = f(4, time_points);
	mulbak = f(10, time_points);
	%f_scaled = [cytcc; bakmcl1; abak; dimbak; mulbak; tbimmcl1];	
	f_scaled = [cytcc; bakmcl1; bak; tbimmcl1];	
end

function [f, h1, h2, h3, h4, h5] = run_simulation(model, activator, T_p, t, x)
	[k, ic, names] = load_model46(activator, model);
	%tBid/tBim
	k(1,1) = 10^x(1);
	k(2,1) = 10^x(2);
	%tBid/tBim
	k(3,1) = 10^x(3);
	k(4,1) = 10^x(4);
	k(5,1) = 10^x(5);
	k(6,1) = 10^x(6);
	k(7,1) = 10^x(7);
	k(8,1) = 10^x(8);
	%tBid/tBim
	k(9,1) = 10^x(9);
	k(1,2) = 10^x(10);
	k(2,2) = 10^x(11);
	%tBid/tBim
	k(3,2) = 10^x(12);
	k(4,2) = 10^x(13);
	k(5,2) = 10^x(14);
	k(6,2) = 10^x(15);
	k(7,2) = 10^x(16);
	k(8,2) = 10^x(17);
	%tBid/tBim
	k(9,2) = 10^x(18);
	%ICs
	%Update B_0 and X_0 with prior info, instead of from load_model46
	%ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];
	ic(1) = x(19);
	ic(13) = x(20);
	f = bcl2model(k, ic, names, t, T_p, activator);
end
