function [k, ic, names, D, u] = load_model(activator, model)
	% load_model		Return kinetic parameters, initial conditions and names of species modelled
	%					for a number of different models. These include:
	%					-	reversible_exp4: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Includes Bak multimers.
	%					-	reversible_exp4_diff: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Diffusion-limited parameter set small enough to be noticeable.
	%						Includes Bak multimers.
	%					-	reversible_exp4_BIAcore: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from BIAcore data and intermediate defaults.
	%						All complexes are reversible. Includes Bak multimers.		
	%					-	reversible_exp4_BIAcore_diff: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from BIAcore data and intermediate defaults.
	%						All complexes are reversible. Diffusion-limited parameter set small enough
	%						to be noticeable. Includes Bak multimers.
	%					-	reversible_exp4_production: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Includes Bak multimers and production and degradation vectors.
	%					-	reversible_exp4_indirect: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Includes Bak multimers. Model without direct activation and with
	%						Bcl-xl.
	%					-	reversible_exp4_diff_indirect: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Diffusion-limited parameter set small enough to be noticeable.
	%						Includes Bak multimers. Model without direct activation and with
	%						Bcl-xl.
	%					-	reversible_exp4_BIAcore_indirect: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from BIAcore data and intermediate defaults.
	%						All complexes are reversible. Includes Bak multimers. Model without direct activation and with
	%						Bcl-xl.		
	%					-	reversible_exp4_BIAcore_diff_indirect: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from BIAcore data and intermediate defaults.
	%						All complexes are reversible. Diffusion-limited parameter set small enough
	%						to be noticeable. Includes Bak multimers. Model without direct activation and with
	%						Bcl-xl.
	%					-	reversible_exp4_production: Initial conditions for modelling experiment 4.
	%						Kinetic parameters estimates from densitometry data. All complexes are
	%						reversible. Includes Bak multimers and production and degradation vectors.
	%						Model without direct activation and with Bcl-xl.
	%
	% Usage:
	%				[k, ic, names, D] = load_model(activator, model)
	%				[k, ic, names, p, u] = load_model(activator, model)
	%
	% Input:
	%				activator = one of 'tBid' or 'tBim'
	%				model = one of 'reversible_exp4', 'reversible_exp4_diff', 'reversible_exp4_BIAcore',
	%						'reversible_exp4_BIAcore_diff'
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
	
	if strcmp(model, 'reversible_exp4')
		[k, ic] = reversible_exp4(activator);
	elseif strcmp(model, 'reversible_exp4_diff')
		[k, ic, D] = reversible_exp4_diff(activator);
	elseif strcmp(model, 'reversible_exp4_BIAcore')
		[k, ic] = reversible_exp4_BIAcore(activator);
	elseif strcmp(model, 'reversible_exp4_BIAcore_diff')
		[k, ic, D] = reversible_exp4_BIAcore_diff(activator);	
	elseif strcmp(model, 'reversible_exp4_production')
		[k, ic, p, u] = reversible_exp4_production(activator);
	elseif strcmp(model, 'reversible_exp4_indirect')
		[k, ic] = reversible_exp4_indirect(activator);
	elseif strcmp(model, 'reversible_exp4_diff_indirect')
		[k, ic, D] = reversible_exp4_diff_indirect(activator);
	elseif strcmp(model, 'reversible_exp4_BIAcore_indirect')
		[k, ic] = reversible_exp4_BIAcore_indirect(activator);
	elseif strcmp(model, 'reversible_exp4_BIAcore_diff_indirect')
		[k, ic, D] = reversible_exp4_BIAcore_diff_indirect(activator);	
	elseif strcmp(model, 'reversible_exp4_production_indirect')
		[k, ic, D, u] = reversible_exp4_production_indirect(activator);
	elseif strcmp(model, 'reversible_exp4_both')
		[k, ic] = reversible_exp4_both(activator);
	else
		throw(MException('ModelError:invalidModel', 'Invalid model: see help'));
	end
end

function [k, ic] = reversible_exp4(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k1a = 10^(-6.00);		 		%Bak activation (1/nM/s)
		k3a = 10^(-4.58);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-4.78);			    %tBim:Mcl-1 dissoc (1/s)	
	else
		k1a = 10^(-3.52);		 		%Bak activation (1/nM/s)
		k3a = 10^(-4.18);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.38);			    %tBid:Mcl-1 dissoc (1/s)
	end
	k1d = 10^(-6);  					%Bak deactivation (1/s)
	k2a = 10^(-3.19);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.19);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-1.63);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-1.28);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-3.13);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-2.79);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-2.02);					%Bak multimer dissociation rate (1/s)
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

function [k, ic] = reversible_exp4_BIAcore(activator)
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

function [k, ic, D] = reversible_exp4_diff(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k1a = 10^(-6.00);		 		%Bak activation (1/nM/s)
		k3a = 10^(-4.58);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-4.78);			    %tBim:Mcl-1 dissoc (1/s)	
	else
		k1a = 10^(-3.52);		 		%Bak activation (1/nM/s)
		k3a = 10^(-4.18);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.38);			    %tBid:Mcl-1 dissoc (1/s)
	end
	k1d = 10^(-6);  					%Bak deactivation (1/s)
	k2a = 10^(-3.19);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.19);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-1.63);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-1.28);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-3.13);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-2.79);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-2.02);					%Bak multimer dissociation rate (1/s)
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

function [k, ic, D] = reversible_exp4_BIAcore_diff(activator)
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
	
	D = 10^(-9);
end

function [k, ic, p, u] = reversible_exp4_production(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k1a = 10^(-6.00);		 		%Bak activation (1/nM/s)
		k3a = 10^(-4.58);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-4.78);			    %tBim:Mcl-1 dissoc (1/s)	
	else
		k1a = 10^(-3.52);		 		%Bak activation (1/nM/s)
		k3a = 10^(-4.18);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.38);			    %tBid:Mcl-1 dissoc (1/s)
	end
	k1d = 10^(-6);  					%Bak deactivation (1/s)
	k2a = 10^(-3.19);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-3.19);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-1.63);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-1.28);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-3.13);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-2.79);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-2.02);					%Bak multimer dissociation rate (1/s)
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

function [k, ic] = reversible_exp4_indirect(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k3a = 10^(-1.16);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-1.05);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-2.30);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-3.02);				%tBim:Bcl-xl dissoc (1/s)
	else
		k3a = 10^(-0.76);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(1.350);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-2.70);				%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-2.11);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1a = 0;							%Bak activation (1/nM/s)
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-2.75);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-2.52);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-3.10);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-2.77);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-3.13);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-4.27);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-4.41);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-1.54);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-3.11);					%Bak:Bcl-xl dissociation (1/s)
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

function [k, ic] = reversible_exp4_BIAcore_indirect(activator)
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

function [k, ic] = reversible_exp4_diff_indirect(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k3a = 10^(-1.16);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-1.05);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-2.30);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-3.02);				%tBim:Bcl-xl dissoc (1/s)
	else
		k3a = 10^(-0.76);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(1.350);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-2.70);				%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-2.11);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1a = 0;							%Bak activation (1/nM/s)
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-2.75);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-2.52);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-3.10);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-2.77);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-3.13);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-4.27);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-4.41);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-1.54);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-3.11);					%Bak:Bcl-xl dissociation (1/s)
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
	k1a = 0;				 		%Bak activation (1/nM/s)
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
		k3a = 10^(-1.16);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-1.05);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-2.30);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-3.02);				%tBim:Bcl-xl dissoc (1/s)
	else
		k3a = 10^(-0.76);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(1.350);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-2.70);				%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-2.11);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1a = 0;							%Bak activation (1/nM/s)
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-2.75);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-2.52);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-3.10);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-2.77);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-3.13);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-4.27);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-4.41);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-1.54);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-3.11);					%Bak:Bcl-xl dissociation (1/s)
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

function [k, ic] = reversible_exp4_both(activator)
	%Kinetic parameters
	if strcmp(activator,'tBim')
		k1a = 10^(-6);					%Bak activation (1/nM/s)
		k3a = 10^(-4.42);			    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-4.61);			    %tBim:Mcl-1 dissoc (1/s)	
		k9a = 10^(-5.23);				%tBim:Bcl-xl assoc (1/nM/s) 
		k9d = 10^(-2.15);				%tBim:Bcl-xl dissoc (1/s)
	else
		k1a = 10^(-3.89);				%Bak activation (1/nM/s)
		k3a = 10^(-4.02);			    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 10^(-2.21);			    %tBid:Mcl-1 dissoc (1/s)
		k9a = 10^(-6);					%tBid:Bcl-xl assoc (1/nM/s)
		k9d = 10^(-1.24);				%tBid:Bcl-xl dissoc (1/s)
	end
	k1d = 10^(-5);  					%Bak deactivation (1/s)
	k2a = 10^(-2.88);  					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 10^(-2.97);  					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-3.33);					%Bak dimerisation (1/nM/s)
	k4d = 10^(-4.79);					%Bak dimerisation dissociation (1/s)
	k5a = 10^(-3.22);					%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-3.44);					%Bak multimer association rate (1/nM/s)
	k6d = 10^(-4.44);					%Bak multimer dissociation rate (1/s)
	k7a = 10^(-3);						%Cytochrome c release via multimer
	k7d = 0;
	k8a = 10^(-1.09);					%Bak:Bcl-xl association (1/nM/s)
	k8d = 10^(-5.99);					%Bak:Bcl-xl dissociation (1/s)
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
