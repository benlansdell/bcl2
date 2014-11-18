%[fits_a, r2_a, bakmcl1r2_a, tbidmcl1r2_a] = fit_exp6_multi([-4 -4 -4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_tbim', './images/models_lsq_fit/exp6_tbid');
%[fits_b, r2_b, bakmcl1r2_b, tbidmcl1r2_b] = fit_exp6_multi_scaleBIAcore([-4 -4 -4 -4 -4 -4, -4, 0, 0], './images/models_lsq_fit/exp6_scaleBIAcore_tbim', './images/models_lsq_fit/exp6_scaleBIAcore_tbid');
%[fits_c, r2_c, bakmcl1r2_c, tbidmcl1r2_c] = fit_exp6_multi_keepratios([-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_keepratios_tbim', './images/models_lsq_fit/exp6_keepratios_tbid');
%[fits_d, r2_d, bakmcl1r2_d, tbidmcl1r2_d] = fit_exp6_multi_diff([-4 -4 -4 -4 -4 -4 -4 -6], './images/models_lsq_fit/exp6_diff_tbim', './images/models_lsq_fit/exp6_diff_tbid');
%[fits_e, r2_e, bakmcl1r2_e, tbidmcl1r2_e] = fit_exp6_multi_indirect([-4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_indirect_tbim', './images/models_lsq_fit/exp6_indirect_tbid');
%[fits_f, r2_f, bakmcl1r2_f, tbidmcl1r2_f] = fit_exp6_multi_scaleBIAcore_indirect([-4 -4 -4 -4, -4, 0, 0], './images/models_lsq_fit/exp6_scaleBIAcore_indirect_tbim', './images/models_lsq_fit/exp6_scaleBIAcore_indirect_tbid');
%[fits_g, r2_g, bakmcl1r2_g, tbidmcl1r2_g] = fit_exp6_multi_keepratios_indirect([-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_keepratios_indirect_tbim', './images/models_lsq_fit/exp6_keepratios_indirect_tbid');

%[fits_h, r2_h, bakmcl1r2_h, tbidmcl1r2_h] = fit_exp6_multi_diff_indirect([-4 -4 -4 -4 -4 -6], './images/models_lsq_fit/exp6_diff_indirect_tbim', './images/models_lsq_fit/exp6_diff_indirect_tbid');
%[fits_i, r2_i, bakmcl1r2_i, tbidmcl1r2_i] = fit_exp6_multi_scaleBIAcore_diff([-4 -4 -4 -4 -4 -4, -4, 0, 0, -6], './images/models_lsq_fit/exp6_scaleBIAcore_diff_tbim', './images/models_lsq_fit/exp6_scaleBIAcore_diff_tbid');
%[fits_j, r2_j, bakmcl1r2_j, tbidmcl1r2_j] = fit_exp6_multi_keepratios_diff([-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -6], './images/models_lsq_fit/exp6_keepratios_diff_tbim', './images/models_lsq_fit/exp6_keepratios_diff_tbid');
%[fits_k, r2_k, bakmcl1r2_k, tbidmcl1r2_k] = fit_exp6_multi_scaleBIAcore_diff_indirect([-4 -4 -4 -4, -4, 0, 0, -6], './images/models_lsq_fit/exp6_scaleBIAcore_diff_indirect_tbim', './images/models_lsq_fit/exp6_scaleBIAcore_diff_indirect_tbid');
%[fits_l, r2_l, bakmcl1r2_l, tbidmcl1r2_l] = fit_exp6_multi_keepratios_diff_indirect([-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -6], './images/models_lsq_fit/exp6_keepratios_diff_indirect_tbim', './images/models_lsq_fit/exp6_keepratios_diff_indirect_tbid');

[fits_m, r2_m, bakmcl1r2_m, tbidmcl1r2_m] = fit_exp6_multi_keepratios_both([-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_keepratios_both_tbim', './images/models_lsq_fit/exp6_keepratios_both_tbid')
[fits_n, r2_n, bakmcl1r2_n, tbidmcl1r2_n] = fit_exp6_multi_both([-4 -4 -4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_both_tbim', './images/models_lsq_fit/exp6_both_tbid')
[fits_o, r2_o, bakmcl1r2_o, tbidmcl1r2_o] = fit_exp6_multi_scaleBIAcore_both([-4 -4 -4 -4 -4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_scaleBIAcore_both_tbim', './images/models_lsq_fit/exp6_scaleBIAcore_both_tbid')
[fits_p, r2_p, bakmcl1r2_p, tbidmcl1r2_p] = fit_exp6_multi_keepratios_diff_both([-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -6], './images/models_lsq_fit/exp6_keepratios_diff_both_tbim', './images/models_lsq_fit/exp6_keepratios_diff_both_tbid')
[fits_q, r2_q, bakmcl1r2_q, tbidmcl1r2_q] = fit_exp6_multi_diff_both([-4 -4 -4 -4 -4 -4 -4 -6], './images/models_lsq_fit/exp6_diff_both_tbim', './images/models_lsq_fit/exp6_diff_both_tbid')
[fits_r, r2_r, bakmcl1r2_r, tbidmcl1r2_r] = fit_exp6_multi_scaleBIAcore_diff_both([-4 -4 -4 -4 -4 -4 -4 -4 -4 -6], './images/models_lsq_fit/exp6_scaleBIAcore_diff_both_tbim', './images/models_lsq_fit/exp6_scaleBIAcore_both_tbid')
[fits_s, r2_s, bakmcl1r2_s, tbidmcl1r2_s] = fit_exp6_multi_keepratios_indirect([-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_keepratios_indirect_13Bak_tbim', './images/models_lsq_fit/exp6_keepratios_indirect_13Bak_tbid')
[fits_t, r2_t, bakmcl1r2_t, tbidmcl1r2_t] = fit_exp6_multi_keepratios_both([-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_keepratios_both_13Bak_tbim', './images/models_lsq_fit/exp6_keepratios_both_13Bak_tbid')
[fits_u, r2_u, bakmcl1r2_u, tbidmcl1r2_u] = fit_exp6_multi_indirect([-4 -4 -4 -4 -4], './images/models_lsq_fit/exp6_indirect_13Bak_tbim', './images/models_lsq_fit/exp6_indirect_13Bak_tbid')
