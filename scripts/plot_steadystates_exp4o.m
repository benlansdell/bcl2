function plot_steadystates_exp4o(filename_tbim, filename_tbid)
	% plot_steadystates_exp4o	Plot equilibrium concentrations as a function of total tBim/tBid spiked
	%			as Experiment 4 proceeds. 
	%
	% Usage:
	%			plot_steadystates_exp4o(filename_tbim, filename_tbid)
	%
	% Input:
	%			filename_tbim = image filename for tbim simulation
	%			filename_tbid = image filename for tbid simluation
	%
	% Example(s):
	%			plot_steadystates_exp4o('./images/test-tbim.eps', './images/test-tbid.eps');

	experiment = 'reversible_exp4';
	T0 = 1:300;
	disp('Computing steady-states for tBid');
	tBid = computeSteadyStateMatrix(experiment, 'tBid', T0);
	disp('Computing steady-states for tBid');
	tBIM = computeSteadyStateMatrix(experiment, 'tBim', T0);

	plotsteadystate(experiment, 'tBim', tBIM, filename_tbim);
	plotsteadystate(experiment, 'tBid', tBid, filename_tbid);
end

function ssdata = computeSteadyStateMatrix(model, activator, T0list)
	ssdata = [];
	[k, ic] = load_model(activator, model);
	for T0 = T0list
		ic(2) = T0;
		ss = steadystate_conc(k, ic)
		ssdata = [ssdata, ss'];
	end
end

function plotsteadystate(model, activator, ssdata, filename)

	len = length(ssdata);
	time_points = [0 10 20 30 40 50 60 70 80];
	B0 = 10;
	M0 = 20;
	window = 60;
	%Plot fit
	[k, ic, names] = load_model(activator, model);
	t = 0:(120*60);
	T_p = T_p_exp4(activator);
	[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p, activator);

	T_added_ = @(t) 0.3 + 0.7*min(max(0, (t-600+window)/window), 1) + 2*min(max(0, (t-1200+window)/window), 1) + 7*min(max(0, (t-1800+window)/window), 1) + 20*min(max(0, (t-2400+window)/window), 1) + 70*min(max(0, (t-3000+window)/window), 1);
	%Additional tBid added to the experiment at 60 minutes...
	if strcmp(activator,'tBid')
		F = @(t) T_added_(t) + 200*min(max(0, (t-3600+window)/window), 1);
	else
		F = @(t) T_added_(t);
	end
	
	%Plot points which were fitted to
	cytcCidx = cellidx(names, 'CytCc');
	cytcMidx = cellidx(names, 'CytCm');
	bakmcl1idx = cellidx(names, 'aBakMcl1');
	tbimmcl1idx = cellidx(names, 'tBimMcl1');
	tbimidx = cellidx(names, 'tBim');
	bakidx = cellidx(names, 'Bak');
	mcl1idx = cellidx(names, 'Mcl1');
	abakidx = cellidx(names, 'aBak');
	dimbakidx = cellidx(names, 'DimBak');

	bakmcl1 = f(bakmcl1idx,:);
	abak = f(abakidx,:);
	bak = f(bakidx,:);
	mcl1 = f(mcl1idx,:);
	tbimmcl1 = f(tbimmcl1idx,:);
	tbim = f(tbimidx,:);
	dimbak = f(dimbakidx,:);
	cytcC = f(cytcCidx,:);
	cytcM = f(cytcMidx,:);
	%b_ss, t_ss, m_ss, d_ss, b_hat_ss, t_hat_ss, m_hat_ss
	colororder = get(0, 'DefaultAxesColorOrder');
	%[b, t, m, d, c, bh, th, mh, ch, b4];
	figure(h1);
	line(t/60, ssdata(1,ceil(F(t))), 'Color', colororder(1,:), 'LineStyle', '--'); %Bak
	line(t/60, ssdata(6,ceil(F(t))), 'Color', colororder(2,:), 'LineStyle', '--'); %aBak
	line(t/60, ssdata(8,ceil(F(t))), 'Color', colororder(3,:), 'LineStyle', '--'); %Mhat
	line(t/60, ssdata(4,ceil(F(t))), 'Color', colororder(4,:), 'LineStyle', '--'); %DimBak
	line(t/60, ssdata(10,ceil(F(t))), 'Color', colororder(5,:), 'LineStyle', '--'); %MultiBak
	legend('Bak', 'Bak*', 'Bak:Mcl-1', 'Dim. Bak*', 'Multi. Bak*', 'Location', 'NorthEast');
	ylim([0 12]);
	%legend('Bak', 'Location', 'NorthEast');
	xlabel('time (min)');
	ylabel('concentration (nM)');
	saveplot(gcf, [filename '-bak.eps']);

	figure(h2);
	line(t/60, ssdata(2,ceil(F(t))), 'Color', colororder(6,:), 'LineStyle', '--') %tBim
	line(t/60, ssdata(7,ceil(F(t))), 'Color', colororder(7,:), 'LineStyle', '--') %That
	ylim([0 23]);
	xlabel('time (min)');
	ylabel('concentration (nM)');
	if (strcmp(activator, 'tBid'))
		legend(activator, [activator ':Mcl-1'], 'Location', 'NorthEast');
	else
		legend(activator, [activator ':Mcl-1'], 'Location', 'East');
	end
	saveplot(gcf, [filename '-activator.eps']);

	figure(h3);
	line(t/60, ssdata(3,ceil(F(t))), 'Color', colororder(8,:), 'LineStyle', '--'); %Mcl-1
	line(t/60, ssdata(8,ceil(F(t))), 'Color', colororder(3,:), 'LineStyle', '--'); %Mhat
	line(t/60, ssdata(7,ceil(F(t))), 'Color', colororder(7,:), 'LineStyle', '--'); %That
	ylim([0 23]);
	xlabel('time (min)');
	ylabel('concentration (nM)');
	if (strcmp(activator, 'tBid'))
		legend('Mcl-1', 'Bak*:Mcl-1', [activator ':Mcl-1'], 'Location', 'NorthEast');
	else
		legend('Mcl-1', 'Bak*:Mcl-1', [activator ':Mcl-1'], 'Location', 'East');
	end
	hold off;
	saveplot(gcf, [filename '-mcl-1.eps']);
end