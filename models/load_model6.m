function [k, ic, names, D, u] = load_model6(activator, model)
	% load_model		Return kinetic parameters, initial conditions and names of species modelled
	%					for a number of different models. These include:
	%					-	reversible_exp6: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Includes Bak multimers.
	%					-	reversible_exp6_diff: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Diffusion-limited parameter set small enough to be noticeable.
	%						Includes Bak multimers.
	%					-	reversible_exp6_BIAcore: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from BIAcore data and intermediate defaults.
	%						All complexes are reversible. Includes Bak multimers.		
	%					-	reversible_exp6_BIAcore_diff: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from BIAcore data and intermediate defaults.
	%						All complexes are reversible. Diffusion-limited parameter set small enough
	%						to be noticeable. Includes Bak multimers.
	%					-	reversible_exp6_production: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Includes Bak multimers and production and degradation vectors.
	%					-	reversible_exp6_indirect: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Includes Bak multimers. Model without direct activation and with
	%						Bcl-xl.
	%					-	reversible_exp6_diff_indirect: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Diffusion-limited parameter set small enough to be noticeable.
	%						Includes Bak multimers. Model without direct activation and with
	%						Bcl-xl.
	%					-	reversible_exp6_BIAcore_indirect: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from BIAcore data and intermediate defaults.
	%						All complexes are reversible. Includes Bak multimers. Model without direct activation and with
	%						Bcl-xl.		
	%					-	reversible_exp6_BIAcore_diff_indirect: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from BIAcore data and intermediate defaults.
	%						All complexes are reversible. Diffusion-limited parameter set small enough
	%						to be noticeable. Includes Bak multimers. Model without direct activation and with
	%						Bcl-xl.
	%					-	reversible_exp6_production: Initial conditions for modelling experiment 6.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Includes Bak multimers and production and degradation vectors.
	%						Model without direct activation and with Bcl-xl.
	%					-	random_parameters_both: pick randomly parameters within the range 10^-1 to 10^-6. Set initial conditions
	%						for use in experiment 6.
	%
	% Usage:
	%				[k, ic, names, D] = load_model6(activator, model)
	%				[k, ic, names, p, u] = load_model6(activator, model)
	%
	% Input:
	%				activator = one of 'tBid' or 'tBim'
	%				model = one of 'reversible_exp6', 'reversible_exp6_diff', 'reversible_exp6_BIAcore',
	%						'reversible_exp6_BIAcore_diff'
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
	%				[k, ic, names] = load_model('tBid', 'reversible_exp4');

	if (~strcmp(activator, 'tBim') & ~strcmp(activator, 'tBid'))
		throw(MException('ActError:invalidActivator', 'Invalid activator: activator must be either tBid or tBim'));
	end
	names = {'Bak'; 'tBim'; 'Mcl1'; 'DimBak'; 'CytCm'; 'aBak'; 'tBimMcl1'; 'aBakMcl1'; 'CytCc'; 'MulBak'; 'Bcl-xl'; 'tBimBclxl'; 'aBakBclxl'};
	D = 0;
	u = [];
	
	if strcmp(model, 'reversible_exp6')
		[k, ic] = reversible_exp6(activator);
	elseif strcmp(model, 'reversible_exp6_diff')
		[k, ic, D] = reversible_exp6_diff(activator);
	elseif strcmp(model, 'reversible_exp6_BIAcore')
		[k, ic] = reversible_exp6_BIAcore(activator);
	elseif strcmp(model, 'reversible_exp6_BIAcore_diff')
		[k, ic, D] = reversible_exp4_BIAcore_diff(activator);	
	elseif strcmp(model, 'reversible_exp4_production')
		[k, ic, p, u] = reversible_exp4_production(activator);
	elseif strcmp(model, 'reversible_exp6_indirect')
		[k, ic] = reversible_exp6_indirect(activator);
	elseif strcmp(model, 'reversible_exp6_diff_indirect')
		[k, ic, D] = reversible_exp6_diff_indirect(activator);
	elseif strcmp(model, 'reversible_exp6_BIAcore_indirect')
		[k, ic] = reversible_exp6_BIAcore_indirect(activator);
	elseif strcmp(model, 'reversible_exp6_BIAcore_diff_indirect')
		[k, ic, D] = reversible_exp6_BIAcore_diff_indirect(activator);	
	elseif strcmp(model, 'reversible_exp6_production_indirect')
		[k, ic, D, u] = reversible_exp6_production_indirect(activator);
	elseif strcmp(model, 'reversible_exp6_both')
		[k, ic] = reversible_exp6_both(activator);
	elseif strcmp(model, 'random_parameters_both');
		[k, ic] = random_parameters_both();
	else
		throw(MException('ModelError:invalidModel', 'Invalid model: see help'));
	end
end

function [k, ic] = reversible_exp6(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k1a = 10^(-3.57);		 		%Bak activation (1/nM/s)
		k3a = 10^(-2.02);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.2);			    %tBim:Mcl-1 dissoc (1/s)	
	else
		k1a = 10^(-5.13);		 		%Bak activation (1/nM/s)
		k3a = 10^(-1.62);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(0.20);			    %tBid:Mcl-1 dissoc (1/s)
	end
	k1d = 10^(-6);  					%Bak deactivation (1/s)
	k2a = 10^(-3.1);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.28);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-3.4);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-3.24);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-2.78);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-4.51);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-4.32);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 0;
	k8d = 0;
	k9a = 0;
	k9d = 0;
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;	      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 0];
end

function [k, ic] = reversible_exp6_BIAcore(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k3a = 10^(-2.49);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-3.53);			    %tBim:Mcl-1 dissoc (1/s)	
	else
		k3a = 10^(-2.09);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-1.12);			    %tBid:Mcl-1 dissoc (1/s)
	end
	k1a = 10^(-4);				 		%Bak activation (1/nM/s)
	k1d = 10^(-6);  					%Bak deactivation (1/s)
	k2a = 10^(-2.86);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.10);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-4.00);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-4.00);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-4.00);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-2.00);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-2.00);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 0;
	k8d = 0;
	k9a = 0;
	k9d = 0;
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 0];
end

function [k, ic, D] = reversible_exp6_diff(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k1a = 10^(-3.57);		 		%Bak activation (1/nM/s)
		k3a = 10^(-2.02);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.2);			    %tBim:Mcl-1 dissoc (1/s)	
	else
		k1a = 10^(-5.13);		 		%Bak activation (1/nM/s)
		k3a = 10^(-1.62);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(0.20);			    %tBid:Mcl-1 dissoc (1/s)
	end
	k1d = 10^(-6);  					%Bak deactivation (1/s)
	k2a = 10^(-3.1);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.28);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-3.4);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-3.24);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-2.78);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-4.51);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-4.32);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 0;
	k8d = 0;
	k9a = 0;
	k9d = 0;
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 0];
	
	D = 10^(-9);
end

function [k, ic, D] = reversible_exp6_BIAcore_diff(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k3a = 10^(-2.49);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-3.53);			    %tBim:Mcl-1 dissoc (1/s)	
	else
		k3a = 10^(-2.09);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-1.12);			    %tBid:Mcl-1 dissoc (1/s)
	end
	k1a = 10^(-4);				 		%Bak activation (1/nM/s)
	k1d = 10^(-6);  					%Bak deactivation (1/s)
	k2a = 10^(-2.86);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.10);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-4.00);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-4.00);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-4.00);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-2.00);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-2.00);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 0;
	k8d = 0;
	k9a = 0;
	k9d = 0;
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 6
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 0];
	
	D = 10^(-9);
end

function [k, ic, p, u] = reversible_exp6_production(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k1a = 10^(-3.57);		 		%Bak activation (1/nM/s)
		k3a = 10^(-2.02);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.2);			    %tBim:Mcl-1 dissoc (1/s)	
	else
		k1a = 10^(-5.13);		 		%Bak activation (1/nM/s)
		k3a = 10^(-1.62);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(0.20);			    %tBid:Mcl-1 dissoc (1/s)
	end
	k1d = 10^(-6);  					%Bak deactivation (1/s)
	k2a = 10^(-3.1);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.28);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-3.4);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-3.24);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-2.78);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-4.51);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-4.32);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 0;
	k8d = 0;
	k9a = 0;
	k9d = 0;
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 0];
	
	%Production rates
	p1 = 10^(-4);
	p = [p1 p1 p1 0 0 0 0 0 0 0 0 0 0];
	%Degradation rates
	u1 = 10^(-5);
	u = [u1 u1 u1 u1 u1 u1 u1 u1 0 0 u1 u1 u1];
end

function [k, ic] = reversible_exp6_indirect(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k3a = 10^(-2.23);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.64);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-4.97);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-2.97);				%tBim:Bcl-xl dissoc (1/s)
	else
		k3a = 10^(-1.83);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(0.24);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-5.37);				%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-2.06);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1a = 0;							%Bak activation (1/nM/s)
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-2.02);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-2.07);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-1.3);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-1.02);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-2.7);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-3.1);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-1.69);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-1.04);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-4.65);					%Bak:Bcl-xl dissociation (1/s)
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
%	B_0 = 7;						%Initial Bak concentration (nM)
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	X_0 = 3;						%Initial Bak:Bcl-xl concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];
end

function [k, ic] = reversible_exp6_BIAcore_indirect(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k3a = 10^(-2.49);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-3.53);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-2.41);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-3.32);				%tBim:Bcl-xl dissoc (1/s)
	else
		k3a = 10^(-2.09);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-1.12);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-2.81);				%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-2.41);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1a = 0;					 		%Bak activation (1/nM/s)
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-2.86);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.10);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-4.00);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-4.00);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-4.00);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-2.00);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-2.00);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-2.56);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-1.99);					%Bak:Bcl-xl dissociation (1/s)
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
%	B_0 = 7;						%Initial Bak concentration (nM)
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	X_0 = 3;						%Initial Bak:Bcl-xl concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];
end

function [k, ic] = reversible_exp6_diff_indirect(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k3a = 10^(-2.23);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.64);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-4.97);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-2.97);				%tBim:Bcl-xl dissoc (1/s)
	else
		k3a = 10^(-1.83);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(0.24);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-5.37);				%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-2.06);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1a = 0;							%Bak activation (1/nM/s)
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-2.02);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-2.07);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-1.3);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-1.02);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-2.7);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-3.1);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-1.69);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-1.04);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-4.65);					%Bak:Bcl-xl dissociation (1/s)
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
%	B_0 = 7;						%Initial Bak concentration (nM)
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	X_0 = 3;						%Initial Bak:Bcl-xl concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];

	D = 10^(-9);
end

function [k, ic, D] = reversible_exp4_BIAcore_diff_indirect(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k3a = 10^(-2.49);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-3.53);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-2.41);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-3.32);				%tBim:Bcl-xl dissoc (1/s)
	else
		k3a = 10^(-2.09);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-1.12);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-2.81);				%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-2.41);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1a = 0;				 			%Bak activation (1/nM/s)
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-2.86);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.10);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-4.00);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-4.00);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-4.00);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-2.00);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-2.00);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-2.56);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-1.99);					%Bak:Bcl-xl dissociation (1/s)
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
%	B_0 = 7;						%Initial Bak concentration (nM)
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	X_0 = 3;						%Initial Bak:Bcl-xl concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];
	
	D = 10^(-9);
end

function [k, ic, p, u] = reversible_exp4_production_indirect(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k3a = 10^(-2.23);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.64);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-4.97);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-2.97);				%tBim:Bcl-xl dissoc (1/s)
	else
		k3a = 10^(-1.83);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(0.24);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-5.37);				%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-2.06);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1a = 0;							%Bak activation (1/nM/s)
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-2.02);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-2.07);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-1.3);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-1.02);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-2.7);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-3.1);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-1.69);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-1.04);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-4.65);					%Bak:Bcl-xl dissociation (1/s)
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
%	B_0 = 7;						%Initial Bak concentration (nM)
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	X_0 = 3;						%Initial Bak:Bcl-xl concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];
	
	%Production rates
	p1 = 10^(-4);
	p = [p1 p1 p1 0 0 0 0 0 0 0 0 0 0];
	%Degradation rates
	u1 = 10^(-5);
	u = [u1 u1 u1 u1 u1 u1 u1 u1 0 0 u1 u1 u1];
end

function [k, ic] = reversible_exp6_both(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k1a = 10^(-2.95);					%Bak activation (1/nM/s)
		k3a = 10^(-2.64);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-3.00);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-1.84);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-3.16);				%tBim:Bcl-xl dissoc (1/s)
	else
		k1a = 10^(-3.91);				%Bak activation (1/nM/s)
		k3a = 10^(-2.24);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-0.60);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-2.24);				%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-2.25);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-3.12);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.42);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-3.38);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-3.18);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-5.60);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-4.09);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-2.92);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-3.16);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-6);					%Bak:Bcl-xl dissociation (1/s)
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

	%Initial conditions for experiment 4
%	B_0 = 7;						%Initial Bak concentration (nM)
	B_0 = 10;						%Initial Bak concentration (nM)
	T_0 = 0.3;						%Initial tBim/tBid concentration (nM)
	M_0 = 20;      					%Initial Mcl-1 concentation (nM)
	C_0 = 5;						%Initial cyt c_m concentration (nM)
	X_0 = 3;						%Initial Bak:Bcl-xl concentration (nM)
	ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];
end

function [k, ic] = random_parameters_both()
	seeds = 10.^(5*rand(16,1)-6);
        %Kinetic parameters
        k1a = seeds(1);                                       %Bak activation (1/nM/s)
        k3a = seeds(2);                           %tBim:Mcl-1 assoc (1/nM/s)
        k3d = seeds(3);                           %tBim:Mcl-1 dissoc (1/s)
        k9a = seeds(4);                               %tBim:Bcl-xl assoc (1/nM/s) 
        k9d = seeds(5);                               %tBim:Bcl-xl dissoc (1/s)
        k1d = seeds(6);                                          %Bak deactivation (1/s)
        k2a = seeds(7);                                       %Bak:Mcl-1 assoc (1/nM/s)
        k2d = seeds(8);                                       %Bak:Mcl-1 dissoc (1/s)
        k4a = seeds(9);                                       %Bak dimerisation (1/nM/s)
        k4d = seeds(10);                                       %Bak dimerisation dissociation (1/s)
        k5a = seeds(11);                                       %Bak auto-activation (1/nM/s)
        k5d = 0;
        k6a = seeds(12);                                       %Bak multimer association rate (1/nM/s)
        k6d = seeds(13);                                       %Bak multimer dissociation rate (1/s)
        k7a = seeds(14);                                          %Cytochrome c release via multimer
        k7d = 0;
        k8a = seeds(15);                                       %Bak:Bcl-xl association (1/nM/s)
        k8d = seeds(16);                                  %Bak:Bcl-xl dissociation (1/s)
        k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d; k7a k7d; k8a k8d; k9a k9d];

        %Initial conditions for experiment 4
        B_0 = 10;                                               %Initial Bak concentration (nM)
        T_0 = 0.3;                                              %Initial tBim/tBid concentration (nM)
        M_0 = 20;                                       %Initial Mcl-1 concentation (nM)
        C_0 = 5;                                                %Initial cyt c_m concentration (nM)
        X_0 = 3;                                                %Initial Bak:Bcl-xl concentration (nM)
        ic = [B_0 T_0 M_0 0 C_0 0 0 0 0 0 0 0 X_0];
end
