function T_p = T_p_exp46(activator)
	% T_p_exp46		Return tBid/tBIM spike-in derivative function for modelling Experiment 4. Every 10
	%				minutes increasing levels of protein are spiked-in: 0.3, 1, 3, 10, 30, 100, 300 nM
	%				Levels increase linearly over 60 seconds at each 10 minute interval -- derivative
	%				is a piece-wise constant function.
	%
	% Usage:
	%				T_p = T_p_exp46(activator)
	%
	% Input:
	%				activator = one of 'tBid' or 'tBim'
	%
	% Output:
	%				T_p = function handle for spike-in function. For time, t, in seconds returns
	%						spiked-in concentration of activator at that time. 
	%
	% Examples:
	%				T_p = T_p_exp46('tBid');
	%				T_p(590)
	%					ans =
	%						0.0117

	window = 60;

	%Additional tBid added to the experiment at 60 minutes...
	T_p = @(t) (0.27*(t>600-window).*(t<=600) + 2.7*(t>1200-window).*(t<=1200) + 27*(t>1800-window).*(t<=1800) + 270*(t>2400-window).*(t<=2400))/window;
end
