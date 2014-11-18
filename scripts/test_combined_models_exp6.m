%Script to test the combined activation model without direct and without autoactivation
%to see which is 'more important'...

model = 'reversible_exp6_both';
t = (-50*60):(181*60);

activator = 'tBid';
[k, ic, names] = load_model6(activator, model);
f = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
bakmcl1_tbid = f(8,:);
mcl1_tbid = f(3,:);
tbimmcl1_tbid = f(7,:);
activator = 'tBim';
[k, ic, names] = load_model6(activator, model);
f = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
bakmcl1_tbim = f(8,:);
mcl1_tbim = f(3,:);
tbimmcl1_tbim = f(7,:);
colororder = get(0, 'DefaultAxesColorOrder');

activator = 'tBid';
k(1,1) = 0;
[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
saveplot(h1, './images/demos/test_combined_models_exp6_nodirect-bak_tbid.eps');
saveplot(h2, './images/demos/test_combined_models_exp6_nodirect-act_tbid.eps');
figure(h3);
line(t/60, bakmcl1_tbid, 'Color', colororder(3,:), 'LineStyle', '--');
line(t/60, mcl1_tbid, 'Color', colororder(8,:), 'LineStyle', '--');
line(t/60, tbimmcl1_tbid, 'Color', colororder(7,:), 'LineStyle', '--');
saveplot(h3, './images/demos/test_combined_models_exp6_nodirect-mcl1_tbid.eps');
saveplot(h4, './images/demos/test_combined_models_exp6_nodirect-cytc_tbid.eps');
saveplot(h5, './images/demos/test_combined_models_exp6_nodirect-bclxl_tbid.eps');

activator = 'tBim';
[k, ic, names] = load_model6(activator, model);
k(1,1) = 0;
[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
saveplot(h1, './images/demos/test_combined_models_exp6_nodirect-bak_tbim.eps');
saveplot(h2, './images/demos/test_combined_models_exp6_nodirect-act_tbim.eps');
figure(h3);
line(t/60, bakmcl1_tbim, 'Color', colororder(3,:), 'LineStyle', '--');
line(t/60, mcl1_tbim, 'Color', colororder(8,:), 'LineStyle', '--');
line(t/60, tbimmcl1_tbim, 'Color', colororder(7,:), 'LineStyle', '--');
saveplot(h3, './images/demos/test_combined_models_exp6_nodirect-mcl1_tbim.eps');
saveplot(h4, './images/demos/test_combined_models_exp6_nodirect-cytc_tbim.eps');
saveplot(h5, './images/demos/test_combined_models_exp6_nodirect-bclxl_tbim.eps');

activator = 'tBid';
[k, ic, names] = load_model6(activator, model);
k(5,1) = 0;
[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
saveplot(h1, './images/demos/test_combined_models_exp6_noauto-bak_tbid.eps');
saveplot(h2, './images/demos/test_combined_models_exp6_noauto-act_tbid.eps');
figure(h3);
line(t/60, bakmcl1_tbid, 'Color', colororder(3,:), 'LineStyle', '--');
line(t/60, mcl1_tbid, 'Color', colororder(8,:), 'LineStyle', '--');
line(t/60, tbimmcl1_tbid, 'Color', colororder(7,:), 'LineStyle', '--');
saveplot(h3, './images/demos/test_combined_models_exp6_noauto-mcl1_tbid.eps');
saveplot(h4, './images/demos/test_combined_models_exp6_noauto-cytc_tbid.eps');
saveplot(h5, './images/demos/test_combined_models_exp6_noauto-bclxl_tbid.eps');

activator = 'tBim';
[k, ic, names] = load_model6(activator, model);
k(5,1) = 0;
[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p_exp6(activator), activator);
saveplot(h1, './images/demos/test_combined_models_exp6_noauto-bak_tbim.eps');
saveplot(h2, './images/demos/test_combined_models_exp6_noauto-act_tbim.eps');
figure(h3);
line(t/60, bakmcl1_tbim, 'Color', colororder(3,:), 'LineStyle', '--');
line(t/60, mcl1_tbim, 'Color', colororder(8,:), 'LineStyle', '--');
line(t/60, tbimmcl1_tbim, 'Color', colororder(7,:), 'LineStyle', '--');
saveplot(h3, './images/demos/test_combined_models_exp6_noauto-mcl1_tbim.eps');
saveplot(h4, './images/demos/test_combined_models_exp6_noauto-cytc_tbim.eps');
saveplot(h5, './images/demos/test_combined_models_exp6_noauto-bclxl_tbim.eps');
