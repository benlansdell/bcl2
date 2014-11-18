%tBIM
obs_tbim = [0  0  5 5 5 5 5; %CytC S/N
			10 10 0 0 0 0 0]'; %Mcl-1:Bak

%tBid
obs_tbid = [0 0  0  0  0  0  0; %CytC S/N
			0 10 10 10 10 10 10]'; %Mcl-1:Bak
			
colormap([0 223/255 1; 1 207/255 0; 0.5 0 0]);
bar(1:7, obs_tbim);
ylabel('final concentration (nM)');
xlabel('initial concentration (nM)');
legend('Cyt. c S/N', 'Bak*:Mcl-1', 'location', 'NorthWest');
ylim([0 12]);
set(gca,'XTickLabel',{'0.3','1','3','10','30', '100', '300'})
saveplot(gcf, './images/demos/densitometry_exp3_tbim.eps');
figure;
colormap([0 223/255 1; 1 207/255 0; 0.5 0 0]);
bar(1:7, obs_tbid);
ylim([0 12]);
ylabel('final concentration (nM)');
xlabel('initial concentration (nM)');
legend('Cyt. c S/N', 'Bak*:Mcl-1', 'location', 'NorthWest');
set(gca,'XTickLabel',{'0.3','1','3','10','30', '100', '300'})
saveplot(gcf, './images/demos/densitometry_exp3_tbid.eps');