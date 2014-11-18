%Computes time taken for each method to solve system upto 90 minutes. Do so for 'typical' 
%kinetic parameters and 'extreme' parameters

experiment = direct_exp4g_lsq;
activator = 'tBid';
[k, ic, t, spike_time, window] = experiment.setup_params(activator);
solvers = {@ode45, @ode23, @ode113, @ode15s, @ode23s, @ode23t, @ode23tb};
ttimes = [];
etimes = [];
plot_solution = 0;

%Typical(T)
k(1,1) = 1e-4;
k(4,1) = 1e-4;
k(5,1) = 1e-4;
k(6,1) = 0;
solver = @ode45;
tstart=tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); ttimes = [ttimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode45_typ.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode23;
tstart=tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); ttimes = [ttimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode23_typ.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode113;
tstart=tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); ttimes = [ttimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode113_typ.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode15s;
tstart=tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); ttimes = [ttimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode15s_typ.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode23s;
tstart=tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); ttimes = [ttimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode23s_typ.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode23t;
tstart=tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); ttimes = [ttimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode23t_typ.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode23tb;
tstart=tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); ttimes = [ttimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode23tb_typ.eps'; print(gcf,  filename, '-depsc');
end

%Extreme(E)
k(1,1) = 1e-1;
k(4,1) = 1e-1;
k(5,1) = 1e-1;
k(6,1) = 0;
solver = @ode45;
tstart = tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); etimes = [etimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode45_ext.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode23;
tstart = tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); etimes = [etimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode23_ext.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode113;
tstart = tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); etimes = [etimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode113_ext.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode15s;
tstart = tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); etimes = [etimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode15s_ext.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode23s;
tstart = tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); etimes = [etimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode23s_ext.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode23t;
tstart = tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); etimes = [etimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode23t_ext.eps'; print(gcf,  filename, '-depsc');
end
solver = @ode23tb;
tstart = tic; experiment.direct_exp4(k, ic, t, spike_time, window, activator, plot_solution, 0, solver); etimes = [etimes toc(tstart)];
if (plot_solution)
	set(gcf, 'paperunits', 'inches'); set(gcf, 'papersize', [6 4]); set(gcf, 'paperposition', [0 0 6 4]);
	filename = './images/direct_exp4g_ode23tb_ext.eps'; print(gcf,  filename, '-depsc');
end

solvers
ttimes
etimes