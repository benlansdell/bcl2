%Script to investigate whether bistable behaviour exists in model by varying initial concentrations
%of tBim/tBid spiked into the system for a set of fitted parameters...

%Really needs the 'final' set of fitted parameters to proceed though...
activator = 'tBim';

[k, ic] = load_model(activator, 'reversible_exp4');
k(1,2) = 10^(-3);
k(1,1) = 10^(-8);
%Vary initial tBim/tBid spiked in and view change in steady-states...
T0 = 1:10:300;
k5a = 10.^(-6+(0:20)/5);
steadystates = zeros(length(T0), length(k5a));
for i = 1:length(T0)
	for j = 1:length(k5a)
		i
		j
		k(5,1) = k5a(j);
		ic(2) = T0(i);
		ss = steadystate_conc(k, ic);
		steadystates(i,j) = ss(10);
	end
end

%Plot solutions
surf(log10(k5a), T0, steadystates);
xlabel('k5a');
ylabel('T_0');
colorbar();
saveplot(gcf, './images/bistability_exp4_tbim.eps');
%Stoichiometric matrix
%S = [-1  0  0  0 -1  0  1  0  0  0  0;
%	  0  0 -1  0  0  0  0  0  1  0  0;
%	  0 -1 -1  0  0  0  0  1  1  0  0;
%	  0  0  0  1  0 -2  0  0  0 -1  2;
%	  1 -1  0 -2  1  0 -1  1  0  2  0;
%	  0  0  1  0  0  0  0  0 -1  0  0;
%	  0  1  0  0  0  0  0 -1  0  0  0;
%	  0  0  0  0  0  1  0  0  0  0 -1];
	  
%null(S)

%k(5,1) = 10^(-5.4)
%ic(2) = 1
%ss = steadystate_conc(k, ic)

