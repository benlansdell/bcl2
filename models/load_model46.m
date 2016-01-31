function [k, ic, names, D, u] = load_model14(activator, model)
	% load_model46		Return kinetic parameters, initial conditions and names of species modelled
	%					for a number of different models. These include:
	%					-	reversible_exp46: Initial conditions for modelling experiment 46.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Includes Bak multimers.
	%
	% Usage:
	%				[k, ic, names, D] = load_model46(activator, model)
	%				[k, ic, names, p, u] = load_model46(activator, model)
	%
	% Input:
	%				activator = one of 'tBid' or 'tBim'
	%				model = one of 'reversible_exp14', 'reversible_exp14_diff', 'reversible_exp14_BIAcore',
	%						'reversible_exp14_BIAcore_diff'
	%
	% Output:
	%				k = 9x2 matrix of kinetic parameters depending on model specified
	%				ic = 13x1 vector of initial concentrations depending on model specified
	%				names = 13x1 cell array of names of each species in same order as ic vector
	%				D = diffusion parameter, if applicable. Otherwise return zero. D or p,u or none
	%					are returned, depending on model specified.
	%				p = 13x1 vector of production rates for model with protein synthesis terms
	%				u = 13x1 vector of protein degradation rates for model with degradation
	%
	% Examples:
	%				[k, ic, names] = load_model46('tBid', 'reversible_exp14');

	if (~strcmp(activator, 'tBim') & ~strcmp(activator, 'tBid'))
		throw(MException('ActError:invalidActivator', 'Invalid activator: activator must be either tBid or tBim'));
	end
	names = {'Bak'; 'tBim'; 'Mcl1'; 'DimBak'; 'CytCm'; 'aBak'; 'tBimMcl1'; 'aBakMcl1'; 'CytCc'; 'MulBak'; 'Bcl-xl'; 'tBimBclxl'; 'aBakBclxl'};
	D = 0;
	u = [];
	
	if strcmp(model, 'reversible_exp46')
		[k, ic] = reversible_exp46(activator);
	else
		throw(MException('ModelError:invalidModel', 'Invalid model: see help'));
	end
end

function [k, ic] = reversible_exp46(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k1a = 10^(-2.369);		 		%Bak activation (1/nM/s)
		k3a = 10^(-2.78);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-4.21);			    %tBim:Mcl-1 dissoc (1/s)	
	else
		k1a = 10^(-2.23);		 		%Bak activation (1/nM/s)
		k3a = 10^(-2.37);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-1.81);			    %tBid:Mcl-1 dissoc (1/s)
	end
	k1d = 10^(-6);  					%Bak deactivation (1/s)
	k2a = 10^(-3.34);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.99);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-2.82);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-2.61);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-3.48);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-4.12);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-3.42);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 0;
	k8d = 0;
	k9a = 0;
	k9d = 0;
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 46
	B_0 = 0.89;						%Initial Bak concentration (nM)
	T_0 = 0.03;						%Initial tBim/tBid concentration (nM)
	M_0 = 35;	      				%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	X_0 = 1.15;						%Initial Bcl-xl concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];
end