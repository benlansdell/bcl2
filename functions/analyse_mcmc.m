%Analyse results

%Try different runs, see if any different
nrun = 2;
data_in = ['./gwmcmc-run' num2str(nrun) '.mat'];
load(data_in);

%Try different burn ins.
burnin = 0.3;
nB = floor(burnin*size(models, 3));
nP = size(models, 1);

%Remove burn in 
models = models(:,:,nB:end);
logP = logP(:,:,nB:end);

%Reshape
models = squeeze(reshape(models, 25, [], 1))';

%Plot traces of each parameter chain
post = reshape(sum(logP, 1), 1, []);

%Find max prob
maxidx = find(post == max(post))
post_params = models(maxidx,:);
map_param = post_params(1,:);

%Plot density plots of 2D marginalizations of parameters
for i = 1:nP
	for j = 1:nP
		subplot(nP,nP, (i-1)*nP+j);
		smoothhist2D(models(:,[i,j]), 5, [100, 100], 0.05)
		hold on 
		plot(map_param(i), map_param(j), 'r.')
	end
end
saveplot(gcf, ['./images/pairwise_postdensity_run_' num2str(nrun) '.png'], 'png', [20 20])
%prctile(models',[5 50 95]) 

%Take and simulate the maximum AP... see if it looks any good at all
model = 'reversible_exp46';
tbim_params = [1 3 4 6 7 8 9 10 11 13 14 15 17 18 19 20 21 22 24 25];
tbid_params = [2 3 5 6 7 8 9 10 12 13 14 16 17 18 19 20 21 23 24 25];
fn_out = ['./images/map_sim_tbim_run_' num2str(nrun) '.eps'];
solve_ode('tBim', model, map_param(tbim_params), fn_out);
fn_out = ['./images/map_sim_tbid_run_' num2str(nrun) '.eps'];
solve_ode('tBid', model, map_param(tbid_params), fn_out);

%Simulate a bunch of different samples from the posterior and see what they look like
nS = 100;
for idx = 1:nS
	[f, h1, h2, h3, h4, h5] = solve_ode('tBim', model, map_param(tbim_params), fn_out);
	[f, h1, h2, h3, h4, h5] = solve_ode('tBid', model, map_param(tbid_params), fn_out);

end

%Compare to max prob without a prior and using local optimization method