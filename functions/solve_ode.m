function [f_scaled, h1, h2, h3, h4, h5] = solve_ode(activator, model, x, fn_out)	
	t = (-50*60):(181*60);
	time_points = (50*60) + 60*[0 10 20 30 40 50 60 90 120 150] + 1;
	[f, h1, h2, h3, h4, h5] = run_simulation(model, activator, T_p_exp46(activator), t, x, fn_out);	
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

function [f, h1, h2, h3, h4, h5] = run_simulation(model, activator, T_p, t, x, fn_out)
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
	if length(fn_out) > 0
		[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p, activator);
		%saveplot(h1, './worksheets/2016-02-01_testrun/testsim1.eps')
		%saveplot(h2, './worksheets/2016-02-01_testrun/testsim2.eps')
		%saveplot(h3, './worksheets/2016-02-01_testrun/testsim3.eps')
		%saveplot(h4, './worksheets/2016-02-01_testrun/testsim4.eps')
	else
		f = bcl2model(k, ic, names, t, T_p, activator);	 
	end
end