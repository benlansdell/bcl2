obs_tbim = [ 0.1604616858    0.340797019     0.3178423615    0.646544851     0.6890785763    0.5761470086    0.5300241985    0.6152542908    0.7180002583    1.3706083633    3.3400985393    4.5111404979; %tBim: cytc sn
4.8395383142    4.659202981     4.6821576385    4.353455149     4.3109214237    4.4238529914    4.4699758015    4.3847457092    4.2819997417    3.6293916367    1.6599014607    0.4888595021; %tBim: cytc mitos
0.7130897909    1.2657178698    2.2305084912    5.3911422325    5.9342776664    6.8080349053    6.4185076651    5.5277042771    5.4554974144    3.9428302241    2.1299185181    1.1995677565; %tbim: bak:mcl1
0       0.3     1       3       10      13.1919650947   13.5814923349   14.4722957229   14.5445025856   16.0571697759   17.8700814819   18.8004322435]'; %tBim:Mcl-1
obs_tbid = [0.2659396793    0.4844349888    0.5996917791    0.7916749669    0.8631459308    0.9225380569    1.1589607123    1.1622242212    1.1868011104    1.0337015607    1.4091379022    1.3687224809; %tBid: cytc sn
4.7340603207    4.5155650112    4.4003082209    4.2083250331    4.1368540692    4.0774619431    3.8410392877    3.8377757788    3.8131988896    3.9662984393    3.5908620978    3.6312775191; %tBid: cytc mitos
1.1274725442    3.1473141009    6.0579847758    8.5225353103    7.7688644212    8.9862966602    8.2652913307    8.9411991903    8.353143432     7.7712134342    6.0903681782    5.2575903352; %tBid: bak:Mcl-1
0       0.3     1       3       10      11.0137033398   11.7347086693   11.0588008097   11.646856568    12.2287865658   13.9096318218   14.7424096648]'; %tBid:Mcl-1

clf;
figure;
%colormap([0 223/255 1; 1 207/255 0; 0.5 0 0]);
%bar([0:10:120 150 180], obs_tbim(:,1:2));
%xlabel('time (min)');
%ylabel('concentration (nM)');
%legend('Cyt. c S/N', 'Cyt. c mito', 'location', 'NorthEast');
%ylim([0 7]);
%xlim([-4 184]);
%saveplot(gcf, './images/demos/densitometry_exp14_tbima.eps');

%bar([0:10:120 150 180], obs_tbim(:,3:4));
%xlabel('time (min)');
%ylabel('concentration (nM)');
%legend('Bak*:Mcl-1', 'tBIM:Mcl-1', 'location', 'NorthEast');
%ylim([0 25]);
%xlim([-4 184]);
%saveplot(gcf, './images/demos/densitometry_exp14_tbimb.eps');

%clf;
%figure;
%colormap([0 0 143/255; 0 223/255 1; 1 207/255 0; 0.5 0 0]);
%bar([0:10:120 150 180], obs_tbid(:,1:2));
%xlim([-4 184]);
%ylim([0 6])
%xlabel('time (min)');
%ylabel('concentration (nM)');
%legend('Cyt. c S/N', 'Cyt. c mito', 'location', 'NorthEast');
%saveplot(gcf, './images/demos/densitometry_exp14_tbida.eps');

%figure;
%colormap([0 0 143/255; 0 223/255 1; 1 207/255 0; 0.5 0 0]);
%bar([0:10:120 150 180], obs_tbid(:,3:4));
%xlim([-4 184]);
%ylim([0 16])
%xlabel('time (min)');
%ylabel('concentration (nM)');
%legend('Bak*:Mcl-1', 'tBid:Mcl-1', 'location', 'NorthWest');
%saveplot(gcf, './images/demos/densitometry_exp14_tbidb.eps');


bar([0:10:90 120 150], obs_tbim(:,1:4));
set(gca, 'FontSize', 10)
xlabel('time (min)');
ylabel('concentration (nM)');
legend('Cyt c S/N', 'Cyt c mito', 'Bak*:Mcl-1', 'tBIM:Mcl-1', 'location', 'NorthEast');
ylim([0 30]);
xlim([-4 184]);
saveplot(gcf, './images/demos/densitometry_exp14_tbim_new.eps');

figure;
%colormap([0 0 143/255; 0 223/255 1; 1 207/255 0; 0.5 0 0]);
bar([0:10:90 120 150], obs_tbid(:,1:4));
set(gca, 'FontSize', 10)
xlim([-4 184]);
ylim([0 16])
xlabel('time (min)');
ylabel('concentration (nM)');
legend('Cyt c S/N', 'Cyt c mito', 'Bak*:Mcl-1', 'tBid:Mcl-1', 'location', 'NorthWest');
saveplot(gcf, './images/demos/densitometry_exp14_tbid_new.eps');
