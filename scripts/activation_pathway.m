function [G_d, G_a] = activation_pathway(activator, model)
	% activation_pathway	Function to compute the flux of Bak molecules through both the direct
	%						activation pathway and the auto-activation pathway for a given set of 
	%						kinetic parameters and activating BH3-only protein.
	%
	% Usage:
	%				[G_d, G_a] = activation_pathway(activator, model)
	%
	% Input:
	%				activator = one of 'tBid' or 'tBim'
	%				model = one of set of kinetic parameters, etc, see choices in load_model()
	%
	% Output:
	%				G_d = Evaluation of \int_0^{120*60}k_{1a}B(t)T(t)\,dt. The 'flux' of Bak through
	%					the direct activation pathway for the given model parameters
	%				G_a = Evaluation of \int_0^{120*60}k_{5a}B(t)\hat{B}(t)\,dt. The 'flux' of Bak
	%					through the auto-activation pathway for the given model parameters
	%
	% Examples:
	%				[G_d, G_a] = activation_pathway('tBid', 'reversible_exp4_both');

	[k, ic, names] = load_model(activator, model);
	T_p = T_P_exp4(activator);
	conc = bcl2model(k, ic, names, 0:(120*60), T_p, activator);
	tbim = conc(2,:);
	bak = conc(1,:);
	abak = conc(6,:);
	G_d = sum(k(1,1)*bak.*tbim);
	G_a = sum(k(5,1)*bak.*abak);
end