function nonnegativitytest
	%Script to test the effect enforcing non-negativity has on the time
	%time taken to run lsqnonlin...
	activator = 'tBid';
	%Setup parameters, choose some that converge for both with and without non-negative constraint
	[k, ic, t, spike_time, window] = setup_params(activator);
	k(1,1) = 10^(-2);
	k(4,1) = 10^(-4);
	k(5,1) = 10^(-6);
	%Scale BIAcore association numbers by factor 10^D
	D = 0;
	d = -6;
	k(2,1) = k(2,1)*10^D;
	k(3,1) = k(3,1)*10^D;

	%Without non-negativity constraint
	options = odeset('MaxStep', 5);
	solver = @ode15s;
	tstart = tic;
	for i = 1:10
		direct_exp4(k, ic, 10^d, t, spike_time, window, activator, 1, 0, solver, options);
	end
	t_wo = toc(tstart);

	%With non-negativity constraint
	options = odeset('MaxStep', 5, 'NonNegative', [1 1 1 1 1 1 1 1 1]);
	tstart = tic;
	for i = 1:10
		direct_exp4(k, ic, 10^d, t, spike_time, window, activator, 1, 0, solver, options);
	end
	t_w = toc(tstart);

	disp(['Without non-negative constraints (s): ' num2str(t_wo)]);
	disp(['With non-negative constraints (s): ' num2str(t_w)]);
end

function [k, ic, t, spike_time, window, names] = setup_params(activator)
	if strcmp(activator,'tBim')
		spike_time = 50;
		k3a = 3.20e-3;				    %tBim:Mcl-1 assoc (1/nM/s)
		k3d = 2.96e-4;				    %tBim:Mcl-1 dissoc (1/s)	
	elseif strcmp(activator,'tBid')
		spike_time = 60;
		k3a = 8.10e-3;				    %tBid:Mcl-1 assoc (1/nM/s)
		k3d = 7.52e-2;				    %tBid:Mcl-1 dissoc (1/s)
	else
		throw(MException('ActError:invalidActivator', 'Invalid activator: activator must be either tBid or tBim'));
	end	
	k1a = 10^(-4);				 		%Bak activation (1/nM/s)
	k1d = 0;         					%Bak deactivation (1/s)
	k2a = 1.37e-3;   					%Bak:Mcl-1 assoc (1/nM/s)
	k2d = 7.98e-4;   					%Bak:Mcl-1 dissoc (1/s)
	k4a = 10^(-4);						%Bak dimerisation (1/nM/s)
	k4d = 0;							%Bak dimerisation dissociation (1/s)
	k5a = 10^(-4);						%Bak auto-activation (1/nM/s)
	k5d = 0;
	k6a = 10^(-4);						%Cytochrome c release (1/nM/s)
	k6d = 0;
	k = [k1a k1d; k2a k2d; k3a k3d; k4a k4d; k5a k5d; k6a k6d];
	B_0 = 10;							%Initial Bak concentration (nM)
	T_0 = 0.3;							%Initial tBim/tBid concentration (nM)
	M_0 = 20;        					%Initial Mcl-1 concentation (nM)
	D_0 = 0;							%Initial dim. Bak concentration (nM)
	C_0 = 5;							%Initial cyt c_m concentration (nM)
	ic = [B_0 T_0 M_0 D_0 C_0 0 0 0 0];
	t = 0:5400;      					%Experiment time (s)
	window = 60;						%Not used in experiment 3
	names = {'Bak'; 'tBim'; 'Mcl1'; 'DimBak'; 'CytCm'; 'aBak'; 'tBimMcl1'; 'aBakMcl1'; 'CytCc'};
end

%Diffusion dynamically affects both the association and dissociation rates of solution-membrane
%interactions
%k_f - membrane-solution forward reaction
function k = k_f(k_on, k_off, D, R)
	k_plus = 10^(6)*D;
	k = k_plus*k_on/(k_plus + R*k_on);
end

%k_r - membrane-solution reverse reaction
function k = k_r(k_on, k_off, D, R)
	k_plus = 10^(6)*D;
	k = k_plus*k_off/(k_plus + R*k_on);
end

function conc = direct_exp4 (k, ic, D, t, spike_time, window, activator, plot_result, plot_reduced, solver, options)
	
	if (nargin == 6)
		plot_result = 0;
	end
	
	%New tBim added to system...
	T_added_prime_ = @(t) (0.7*(t>600-window).*(t<=600) + 2*(t>1200-window).*(t<=1200) + 7*(t>1800-window).*(t<=1800) + 20*(t>2400-window).*(t<=2400) + 70*(t>3000-window).*(t<=3000))/window;
	T_added_ = @(t) 0.3 + 0.7*min(max(0, (t-600+window)/window), 1) + 2*min(max(0, (t-1200+window)/window), 1) + 7*min(max(0, (t-1800+window)/window), 1) + 20*min(max(0, (t-2400+window)/window), 1) + 70*min(max(0, (t-3000+window)/window), 1);
	%Additional tBid added to the experiment at 60 minutes...
	if strcmp(activator,'tBid')
		T_added_prime = @(t) T_added_prime_(t) + 200*(t>3600-window).*(t<=3600)/window;
		T_added = @(t) T_added_(t) + 200*min(max(0, (t-3600+window)/window), 1);
	else
		T_added_prime = @(t) T_added_prime_(t);
		T_added = @(t) T_added_(t);
	end
	ic(2) = T_added(0);

	%Activation of Bak
	v1 = @(x, k, ic, t) k_f(k(1,1), k(1,2), D, x(1))*x(1)*x(2);
	%Formation of Bak:Mcl-1 complex
	v2 = @(x, k, ic, t) k_f(k(2,1), k(2,2), D, x(6))*x(3)*x(6) - k_r(k(2,1), k(2,2), D, x(6))*x(8);
	%Formation of tBim(d)/Mcl-1 complex
	v3 = @(x, k, ic, t) k(3,1)*x(2)*x(3) - k(3,2)*x(7);
	%Dimerisation of Bak*
	v4 = @(x, k, ic, t) k(4,1)*x(6)*x(6) - k(4,2)*x(4);
	%Auto-activation of Bak
	v5 = @(x, k, ic, t) k(5,1)*x(6)*x(1);
	%Cytochrome C release
	v6 = @(x, k, ic, t) k(6,1)*x(4)*x(5);
	
	df = @(t, x) [-v1(x, k, ic, t) - v5(x, k, ic, t); -v3(x, k, ic, t) + T_added_prime(t); -v2(x, k, ic, t) - v3(x, k, ic, t); v4(x, k, ic, t); -v6(x, k, ic, t); v1(x, k, ic, t) + v5(x, k, ic, t) - v2(x, k, ic, t) - 2*v4(x, k, ic, t); v3(x, k, ic, t); v2(x, k, ic, t);
					v6(x, k, ic, t)];
	%Needs to be set or may skip over tBid/tBim being added...
	[T, f] = solver(df, t, ic, options);

	len = length(t);
	if (plot_result ~= 0)
		figure
		plot(t/60, f)
		xlabel('time (min)');
		ylabel('concentration (nM)');
		ylim([0 30]);
		title([num2str(ic(1)) 'nM of Bak, increasing levels of ' activator ' and ' num2str(ic(3)) 'nM of Mcl-1']);
		legend('Bak', activator, 'Mcl-1', 'Dim. Bak*', 'Cyt-C_M', 'Bak*', [activator ':Mcl-1'], 'Bak*:Mcl-1', 'Cyt-C_C', -1);
	end

	%Return concentrations of experiment
	conc = f';
end
