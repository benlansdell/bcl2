%Script to test the combined activation model without direct and without autoactivation
%to see which is 'more important'...

model = 'reversible_exp4_both';
t = (-50*60):(120*60);

activator = 'tBid';
[k, ic, names] = load_model(activator, model);
k(1,1) = 0;
[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p_exp4(activator), activator);
saveplot(h1, './images/demos/test_combined_models_nodirect-bak_tbid.eps');
saveplot(h2, './images/demos/test_combined_models_nodirect-act_tbid.eps');
saveplot(h3, './images/demos/test_combined_models_nodirect-mcl1_tbid.eps');
saveplot(h4, './images/demos/test_combined_models_nodirect-cytc_tbid.eps');
saveplot(h5, './images/demos/test_combined_models_nodirect-bclxl_tbid.eps');

activator = 'tBim';
[k, ic, names] = load_model(activator, model);
k(1,1) = 0;
[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p_exp4(activator), activator);
saveplot(h1, './images/demos/test_combined_models_nodirect-bak_tbim.eps');
saveplot(h2, './images/demos/test_combined_models_nodirect-act_tbim.eps');
saveplot(h3, './images/demos/test_combined_models_nodirect-mcl1_tbim.eps');
saveplot(h4, './images/demos/test_combined_models_nodirect-cytc_tbim.eps');
saveplot(h5, './images/demos/test_combined_models_nodirect-bclxl_tbim.eps');

activator = 'tBid';
[k, ic, names] = load_model(activator, model);
k(5,1) = 0;
[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p_exp4(activator), activator);
saveplot(h1, './images/demos/test_combined_models_noauto-bak_tbid.eps');
saveplot(h2, './images/demos/test_combined_models_noauto-act_tbid.eps');
saveplot(h3, './images/demos/test_combined_models_noauto-mcl1_tbid.eps');
saveplot(h4, './images/demos/test_combined_models_noauto-cytc_tbid.eps');
saveplot(h5, './images/demos/test_combined_models_noauto-bclxl_tbid.eps');

activator = 'tBim';
[k, ic, names] = load_model(activator, model);
k(5,1) = 0;
[f, h1, h2, h3, h4, h5] = bcl2model(k, ic, names, t, T_p_exp4(activator), activator);
saveplot(h1, './images/demos/test_combined_models_noauto-bak_tbim.eps');
saveplot(h2, './images/demos/test_combined_models_noauto-act_tbim.eps');
saveplot(h3, './images/demos/test_combined_models_noauto-mcl1_tbim.eps');
saveplot(h4, './images/demos/test_combined_models_noauto-cytc_tbim.eps');
saveplot(h5, './images/demos/test_combined_models_noauto-bclxl_tbim.eps');
