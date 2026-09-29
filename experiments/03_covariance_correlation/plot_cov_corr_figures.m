%% Report figures: posterior covariance and correlation vs noise
%  Uses the results saved by compare_cov_noise.m (run it with
%  useTuned = false and useTuned = true first).
clear; close all; clc;

P = repo_paths();
orig  = load(fullfile(P.results, 'compare_cov_noise_results.mat'));
tuned = load(fullfile(P.results, 'compare_cov_noise_results_tuned.mat'));
noise_levels = tuned.noise_levels;

outFolder = fullfile(P.report, 'cov_corr');
if ~exist(outFolder,'dir') mkdir(outFolder); end

names  = {'width', 'position'};
labels = {'w', 'p'};

for e = 1:2
    E  = tuned.exps(e);
    Eo = orig.exps(e);

    % ---- covariance scaled by sigma^2 (tuned MCMC vs Laplace (F^{-1})_{12}) ----
    fig = figure;
    semilogx(noise_levels, E.cov ./ noise_levels.^2, 'o-', 'LineWidth', 1.5); hold on;
    semilogx(noise_levels([1 end]), E.cov_unit*[1 1], 'k--', 'LineWidth', 1);
    xlabel('$\sigma_{\mathrm{noise}}$', 'Interpreter', 'latex');
    ylabel(['Cov$(a,' labels{e} ')/\sigma^2_{\mathrm{noise}}$'], 'Interpreter', 'latex');
    legend({'MCMC', 'Laplace'}, 'Location', 'southeast');
    title('Scaled posterior covariance', 'Interpreter', 'latex');
    grid on;
    saveTikz(fig, ['cov_vs_noise_' names{e}], outFolder)

    % ---- correlation (original and tuned MCMC vs Laplace) ----
    fig = figure;
    semilogx(noise_levels, Eo.corr, 's:', 'LineWidth', 1.2); hold on;
    semilogx(noise_levels, E.corr, 'o-', 'LineWidth', 1.5);
    semilogx(noise_levels([1 end]), E.corr_lapl*[1 1], 'k--', 'LineWidth', 1);
    ylim([-1 1]);
    xlabel('$\sigma_{\mathrm{noise}}$', 'Interpreter', 'latex');
    ylabel(sprintf('corr$(a,%s)$', labels{e}), 'Interpreter', 'latex');
    legend({'MCMC (original)', 'MCMC (tuned)', 'Laplace'}, 'Location', 'southeast');
    title('Posterior correlation', 'Interpreter', 'latex');
    grid on;
    saveTikz(fig, ['corr_vs_noise_' names{e}], outFolder)
end

% Numbers quoted in the report
for e = 1:2
    fprintf('%s: Laplace corr = %+.3f, Cov/sigma^2 = %.3e, std/sigma = [%.3e %.3e]\n', ...
        names{e}, tuned.exps(e).corr_lapl, tuned.exps(e).cov_unit, tuned.exps(e).std_unit);
end
