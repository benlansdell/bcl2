function T_p = T_p_none(activator)
	% T_p_none		Return tBid/tBIM spike-in derivative function for zero additions of tBid/tBIM.
	%				Just returns the zero function.
	%
	% Usage:
	%				T_p = T_p_exp4(activator)
	%
	% Input:
	%				activator = one of 'tBid' or 'tBim'
	%
	% Output:
	%				T_p = function handle for spike-in function. For time, t, in seconds returns
	%						spiked-in concentration of activator at that time. 
	%
	% Examples:
	%				T_p = T_p_exp4('tBid');
	%				T_p(590)
	%					ans =
	%						0

	T_p = @(t) 0;
end
