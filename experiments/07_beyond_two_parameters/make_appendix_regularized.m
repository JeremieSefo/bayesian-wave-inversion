%% Appendix figure: Gaussian basis (Nb = 80) with flat prior, gradient prior
%  (tau from the evidence) and L2 prior (physical tau = 0.5).
%  Needs regularized_results.mat (from regularized_inversion.m).
clear; close all; clc;

P = repo_paths();
S = load(fullfile(P.data, 'waveData.mat'));
R = load(fullfile(P.results, 'regularized_results.mat'));
outFolder = P.slides;
fs = 14;

n = numel(S.Contrast);  x = linspace(0, 2, n).';
N = numel(S.TrueScattData_noisy);  M = 100;  T = 1.5;
kEff = rel_noise_factor(S.TrueScattData_noisy, T) / sqrt(M);
ct = 0.5 * (S.Contrast + sqrt(S.Contrast.^2 + 4));

jShow = [1 3 6];
rows = {'flat',     'flat prior';
        'gradient', '$\|\rho''\|^2$, $\tau$ from evidence';
        'tikfix',   '$L^2$, physical $\tau = 0.5$'};
fig = figure('Position', [100 100 1400 850]);
tiledlayout(size(rows,1), numel(jShow), 'TileSpacing', 'compact', 'Padding', 'compact');
xm = x <= 0.6;  xx = x(xm).';
for r = 1:size(rows,1)
    for q = 1:numel(jShow)
        j = jShow(q);  sn = R.noise_levels(j);
        Q = R.res.(sprintf('gauss_80_%s_%d', rows{r,1}, j));
        nexttile; hold on;
        fill([xx fliplr(xx)], [Q.lo(xm).' fliplr(Q.hi(xm).')], [0.75 0.8 1], 'EdgeColor', 'none');
        plot(xx, Q.cm(xm), 'b-', 'LineWidth', 1.8);
        plot(xx, ct(xm), 'r--', 'LineWidth', 1.5);
        xlim([0 0.6]);  ylim([0.6 2.0]);
        title(sprintf('%s, $\\sigma = %g$ (%s\\%%)', rows{r,2}, sn, pct_str(100*kEff*sn)), ...
            'Interpreter', 'latex');
        text(0.585, 1.9, sprintf('err %.0f\\%%, coverage %.0f\\%%', 100*Q.err, 100*Q.cov), ...
            'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'FontSize', fs-1);
        if q == 1, ylabel('$c/c_0$', 'Interpreter', 'latex'); end
        if r == size(rows,1), xlabel('$x$', 'Interpreter', 'latex'); end
        if r == 1 && q == 1
            legend({'95\% region', 'posterior mean', 'true'}, 'Interpreter', 'latex', ...
                'Location', 'west', 'FontSize', fs-3);
        end
        grid on; box on; set(gca, 'FontSize', fs);
    end
end
exportgraphics(fig, fullfile(outFolder, 'slide_regularized.png'), 'Resolution', 300);
