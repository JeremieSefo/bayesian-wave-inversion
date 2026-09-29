%% Slide figure: 95% credible-interval width vs noise level (log-log, slopes)
%  CI widths from the production samples of run_slide_results.m
%  (slide_results.mat); Laplace prediction 2*1.96*sigma*std_unit from
%  compare_cov_noise_results_tuned.mat. Widths are relative to the true value.
%  Slope 1 on log-log axes means width proportional to sigma.
clear; close all; clc;

P = repo_paths();
R = load(fullfile(P.results, 'slide_results.mat'));
L = load(fullfile(P.results, 'compare_cov_noise_results_tuned.mat'));
S = load(fullfile(P.data, 'waveData.mat'));
outFolder = P.slides;
fs = 16;
kEff = rel_noise_factor(S.TrueScattData_noisy, 1.5) / sqrt(100);

sig = R.noise_levels;
sfine = logspace(log10(sig(1)), log10(sig(end)), 50);
% (experiment, parameter index, label, true value)
P = {1, 1, 'a', 0.75;  1, 2, 'w', 5e-3;  2, 2, 'p', 7/50};
col = lines(3);

fig = figure('Position', [100 100 1000 560]);
hold on;
h = gobjects(1, size(P,1));  lbl = cell(1, size(P,1));
fprintf('%-3s | %-22s | %-22s | Laplace\n', '', 'slope, all 6 levels', 'slope, sigma >= 0.05');
for k = 1:size(P,1)
    [e, i, name, tv] = P{k,:};
    wdt = arrayfun(@(j) R.R(e,j).q975(i) - R.R(e,j).q025(i), 1:numel(sig)) / tv;
    lap = 2*1.96*sfine*L.exps(e).std_unit(i) / tv;
    c1 = polyfit(log(sig), log(wdt), 1);
    c2 = polyfit(log(sig(2:end)), log(wdt(2:end)), 1);
    fprintf('%-3s | %6.2f                 | %6.2f                 | %.2e\n', name, c1(1), c2(1), ...
        2*1.96*L.exps(e).std_unit(i)/tv);
    loglog(sfine, lap, '--', 'LineWidth', 1.5, 'Color', col(k,:));
    h(k) = loglog(sig, wdt, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', col(k,:), ...
        'MarkerFaceColor', col(k,:));
    if strcmp(name, 'p')
        lbl{k} = sprintf('$%s$: slope %.2f (%.2f for $\\sigma \\ge 0.05$)', name, c1(1), c2(1));
    else
        lbl{k} = sprintf('$%s$: slope %.2f', name, c1(1));
    end
end
set(gca, 'XScale', 'log', 'YScale', 'log');
% position is snapped to the grid: the interval cannot be narrower than one cell
hg = yline((1/3000)/(7/50), ':', 'LineWidth', 1.8, 'Color', col(3,:));
hl = plot(nan, nan, 'k--', 'LineWidth', 1.5);
xlim([0.008 1]);
noise_axis(gca, kEff);
xlabel('noise level $\sigma$ \ (effective relative noise)', 'Interpreter', 'latex');
ylabel('95\% credible-interval width / true value', 'Interpreter', 'latex');
legend([h hl hg], [lbl, {'Laplace: $2\times1.96\,\sigma\sqrt{(G^{-1})_{ii}}$', 'one grid cell (for $p$)'}], ...
    'Interpreter', 'latex', 'Location', 'northwest');
title('Credible-interval width $\propto \sigma$ (slope 1)', 'Interpreter', 'latex');
grid on; box on; set(gca, 'FontSize', fs);

exportgraphics(fig, fullfile(outFolder, 'slide_ciwidth.png'), 'Resolution', 300);
export_tikz_panels(fig, 'slide_ciwidth', fullfile(outFolder, 'tikz'));
