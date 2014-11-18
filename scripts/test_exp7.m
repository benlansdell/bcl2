%Script to run a set of experiment 7 simulations for both direct and indirect activation 
%hypotheses...

t = 0:(120*60);
len = length(t);
T0list = [0.3 1 3 10 30 100 300];
prev_ColorOrder = get(0, 'DefaultAxesColorOrder');
my_ColorOrder = [0.00000   0.00000   0.00000;
			     0.00000   0.00000   0.16000;
			     0.00000   0.00000   0.33000;
  			     0.00000   0.00000   0.50000;
  			     0.00000   0.00000   0.66000;
  				 0.00000   0.00000   0.83000;
  				 0.00000   0.00000   1.00000];
set(0,'DefaultAxesColorOrder',my_ColorOrder);

%Direct model...tbid
cytcc = [];
multibak = [];
bakmcl1 = [];
activator = 'tBid';
[k, ic, names] = load_model6(activator, 'reversible_exp6');
for T0 = T0list
	ic(2) = T0;
	conc = bcl2model(k, ic, names, t, T_p_none, activator);
	cytcc = [cytcc; conc(9,:)];
	multibak = [multibak; conc(10,:)];
	bakmcl1 = [bakmcl1; conc(8,:)];
end
%tBid
obs_tbid = [multibak(:,len)'; %CytC S/N
			bakmcl1(:,len)']; %Mcl-1:Bak
plot(t/60, bakmcl1);
xlabel('time (min)')
ylabel('concentration (nM)');
%colorbar();
saveplot(gcf, './images/demos/test_exp7-direct-tbid.eps');

%Direct model...tbim
cytcc = [];
multibak = [];
bakmcl1 = [];
activator = 'tBim';
[k, ic, names] = load_model6(activator, 'reversible_exp6');
for T0 = T0list
	ic(2) = T0;
	conc = bcl2model(k, ic, names, t, T_p_none, activator);
	cytcc = [cytcc; conc(9,:)];
	multibak = [multibak; conc(10,:)];
	bakmcl1 = [bakmcl1; conc(8,:)];
end
%tBIM
obs_tbim = [multibak(:,len)'; %CytC S/N
			bakmcl1(:,len)']; %Mcl-1:Bak
plot(t/60, bakmcl1);
xlabel('time (min)');
ylabel('concentration (nM)');
%colorbar();
saveplot(gcf, './images/demos/test_exp7-direct-tbim.eps');
			
colormap([0 223/255 1; 1 207/255 0; 0.5 0 0]);
bar(1:7, obs_tbim');
ylabel('final concentration (nM)');
xlabel('initial concentration (nM)');
legend('Multi. Bak', 'Bak*:Mcl-1', 'location', 'NorthWest');
ylim([0 12]);
set(gca,'XTickLabel',{'0.3','1','3','10','30', '100', '300'})
saveplot(gcf, './images/demos/test_exp7_final-direct-tbim.eps');
figure;
colormap([0 223/255 1; 1 207/255 0; 0.5 0 0]);
bar(1:7, obs_tbid');
ylim([0 12]);
ylabel('final concentration (nM)');
xlabel('initial concentration (nM)');
legend('Multi. Bak', 'Bak*:Mcl-1', 'location', 'NorthWest');
set(gca,'XTickLabel',{'0.3','1','3','10','30', '100', '300'})
saveplot(gcf, './images/demos/test_exp7_final-direct-tbid.eps');


%Indirect model...tbid
cytcc = [];
multibak = [];
activator = 'tBid';
[k, ic, names] = load_model6(activator, 'reversible_exp6_indirect');
for T0 = T0list
	ic(2) = T0;
	conc = bcl2model(k, ic, names, t, T_p_none, activator);
	cytcc = [cytcc; conc(9,:)];
	multibak = [multibak; conc(10,:)];
end
%tBid
obs_tbid = [multibak(:,len)'; %CytC S/N
			bakmcl1(:,len)']; %Mcl-1:Bak
plot(t/60, bakmcl1);
xlabel('time (min)')
ylabel('concentration (nM)');
%colorbar();
saveplot(gcf, './images/demos/test_exp7-indirect-tbid.eps');

%Indirect model...tbim
cytcc = [];
multibak = [];
activator = 'tBim';
[k, ic, names] = load_model6(activator, 'reversible_exp6_indirect');
for T0 = T0list
	ic(2) = T0;
	conc = bcl2model(k, ic, names, t, T_p_none, activator);
	cytcc = [cytcc; conc(9,:)];
	multibak = [multibak; conc(10,:)];
end
%tBIM
obs_tbim = [multibak(:,len)'; %CytC S/N
			bakmcl1(:,len)']; %Mcl-1:Bak
plot(t/60, bakmcl1);
xlabel('time (min)')
ylabel('concentration (nM)');
%colorbar();
saveplot(gcf, './images/demos/test_exp7-indirect-tbim.eps');

colormap([0 223/255 1; 1 207/255 0; 0.5 0 0]);
bar(1:7, obs_tbim');
ylabel('final concentration (nM)');
xlabel('initial concentration (nM)');
legend('Multi. Bak', 'Bak*:Mcl-1', 'location', 'NorthWest');
ylim([0 12]);
set(gca,'XTickLabel',{'0.3','1','3','10','30', '100', '300'})
saveplot(gcf, './images/demos/test_exp7_final-indirect-tbim.eps');
figure;
colormap([0 223/255 1; 1 207/255 0; 0.5 0 0]);
bar(1:7, obs_tbid');
ylim([0 12]);
ylabel('final concentration (nM)');
xlabel('initial concentration (nM)');
legend('Multi. Bak', 'Bak*:Mcl-1', 'location', 'NorthWest');
set(gca,'XTickLabel',{'0.3','1','3','10','30', '100', '300'})
saveplot(gcf, './images/demos/test_exp7_final-indirect-tbid.eps');

set(0, 'DefaultAxesColorOrder', prev_ColorOrder);