%Test MCMC code on simple example with only a few parameters

nP = 1000;
sig = 0.2;
prior_sig = 1;
mu = [1, 1];
nW = 300;
nruns = 500;
nPars = 2;
minit = zeros(nPars, nW);

%Write down ODE
%da/dt = kd*c - ka*a*b^2
%db/dt = 2*kd*c - 2*ka*a*b^2
%dc/dt = ka*a*b^2 - kd*c

%True solution 
x0 = [10; 10; 10];
ka = 1;
kd = 2;

%Generate ODE data
df = @(t,x) [kd*x(3)-ka*x(1)*x(2)*x(2); 2*kd*x(3)-2*ka*x(1)*x(2)*x(2); ka*x(1)*x(2)*x(2) - kd*x(3)];
t = linspace(0, .2, nP);
options = odeset('MaxStep', 5, 'NonNegative', [1 1 1 1 1 1 1 1 1 1 1 1 1]);
solver = @ode23t;
[T, f] = solver(df, t, x0, options);

%Add noise 
pts = 1:100:nP;
tpts = T(pts);
data = f(pts,:);
data = data + sig*randn(size(data));

%Plot
plot(T, f, tpts, data, 'o')

%Create likelihood function: SSR between data and simulation with given parameters
ll = @(k) ll_test(k, t, x0, options, data, sig, solver, pts);
%Create prior function 
prior = @(k) prior_test(k, mu, prior_sig);

%Run gwmcmc
logPfuns = {prior ll};	
%make a set of starting points for the entire ensemble of walkers 
for idx = 1:nW
	minit(:,idx) = sampleprior_test(mu, prior_sig); 
end

%Apply the MCMC hammer 
[models,logP]=gwmcmc(minit,logPfuns,nW*nruns); 

%Analyse results
models(:,:,1:floor(size(models,3)*.2))=[]; %remove 20% as burn-in
models=models(:,:)'; %reshape matrix to collapse the ensemble member dimension
scatter(models(:,1),models(:,2))
prctile(models,[5 50 95])
nM = size(models,1);

%Draw some samples from the posterior and plot their solution 
nS = 200;
idxsamp = randsample(1:nM, nS);

samp_models = models(idxsamp, :, 3);
t = linspace(0, .2, nP);
samp_sols = zeros(nS, length(t));

clf
for idx = 1:nS
	hold on 
	ka = samp_models(idx,1);
	kd = samp_models(idx,2);
	df = @(t,x) [kd*x(3)-ka*x(1)*x(2)*x(2); 2*kd*x(3)-2*ka*x(1)*x(2)*x(2); ka*x(1)*x(2)*x(2) - kd*x(3)];
	[T, f] = solver(df, t, x0, options);
	samp_sols(idx, :, :) = f;
	plot(T, f, tpts, data, 'o')
end

%Plot histograms of parameter estimates
