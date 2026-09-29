%% RUN_BOTH_SAMPLERS  Run compare_cov_noise.m with both samplers
%  1. useTuned = false : original fixed-step sampler (as in the report)
%  2. useTuned = true  : proposal scaled to the Laplace posterior width
%  Results: results/compare_cov_noise_results.mat and ..._tuned.mat,
%  used by plot_cov_corr_figures.m and by experiments 04, 05 and 08.
clear;
useTuned = false;  compare_cov_noise;
useTuned = true;   compare_cov_noise;
