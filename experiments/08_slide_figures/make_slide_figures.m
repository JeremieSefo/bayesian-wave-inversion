%% Figures for the Bayesian wave inversion slides (figures/slides)
%  1. slide_model.png     : contrast profiles and their scattered data
%  2. slide_noise.png     : noisy trajectories at two noise levels
%  3. slide_variance.png  : posterior variance ~ sigma^2, correlation flat
%  Needs compare_cov_noise_results_tuned.mat (from compare_cov_noise.m).
clear; close all; clc;
rng(1);

P = repo_paths();
load(fullfile(P.data, 'waveData.mat'));
outFolder = P.slides;
if ~exist(outFolder,'dir') mkdir(outFolder); end
fs = 16;                                   % font size for slides

% Setup as in WAVE_MCMC / WAVEpos_MCMC
nI = 1500 * 2;  n = nI;  Le = 2;  hh = Le / (n-1);
center_freq = 2*pi*0.02e3;
w   = center_freq + 2*0.25*center_freq;
lam = 2*pi / w;
c0 = 1;  r = 2;
separation_mult = 3;  width_scale = 1;
C1_amp = 0.6;  C2_amp = 1;
amp_true = 0.75;  width_true = 5e-3;  pos_true = 7/50;

N  = numel(TrueScattData_noisy);
T  = 1.5;
dt = T/(N-1);
t  = (0:N-1)*dt;
M  = 100;                                  % observed trajectories
kRel = rel_noise_factor(TrueScattData_noisy, T);   % ||sigma W||/||D|| per unit sigma
kEff = kRel / sqrt(M);                     % effective, M trajectories combined
fprintf('relative noise: %.1f*sigma per trajectory, %.2f*sigma effective (M = %d)\n', kRel, kEff, M);

%% ======================
% 1. Contrast profiles and scattered data
% ======================
cases = [amp_true, width_true,   pos_true;
         amp_true, 1.5*width_true, pos_true;
         amp_true, width_true,   pos_true + 0.02];
names = {'true $(a,w,p)$', 'wider $w$', 'shifted $p$'};
styles = {'-', '--', ':'};

fig = figure('Position', [100 100 1200 420]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
for c = 1:size(cases,1)
    [rho, x] = contrast(cases(c,:), separation_mult, width_scale, C1_amp, C2_amp, n, hh, lam, r, c0);
    plot(x, rho, styles{c}, 'LineWidth', 2);
end
xlim([0 0.7]);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('contrast $\rho(x)$', 'Interpreter', 'latex');
legend(names, 'Interpreter', 'latex', 'Location', 'northeast');
title('Medium', 'Interpreter', 'latex');
grid on; box on; set(gca, 'FontSize', fs);

nexttile; hold on;
for c = 1:size(cases,1)
    rho = contrast(cases(c,:), separation_mult, width_scale, C1_amp, C2_amp, n, hh, lam, r, c0);
    plot(t, TrueOper * rho, styles{c}, 'LineWidth', 2);
end
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$D^{\mathrm{sc}}(t) = A\rho$', 'Interpreter', 'latex');
title('Scattered data', 'Interpreter', 'latex');
grid on; box on; set(gca, 'FontSize', fs);

exportgraphics(fig, fullfile(outFolder, 'slide_model.png'), 'Resolution', 300);

%% ======================
% 2. Noisy trajectories X_t = X0 exp(D(t) + sigma W_t)
% ======================
X0 = 1;  nPlot = 10;
sig_show = [5e-3, 5e-2];

fig = figure('Position', [100 100 1200 420]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:2
    dW = sqrt(dt) * randn(nPlot, N-1);
    W  = [zeros(nPlot,1), cumsum(dW, 2)];
    X  = X0 * exp(TrueScattData_noisy.' + sig_show(s)*W);
    nexttile; hold on;
    plot(t, X.', 'LineWidth', 0.8);
    plot(t, X0*exp(TrueScattData_noisy), 'k--', 'LineWidth', 2.5);
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$X_t$', 'Interpreter', 'latex');
    title(sprintf('$\\sigma = %g$: %s\\%% per path, %s\\%% effective', ...
        sig_show(s), pct_str(100*kRel*sig_show(s)), pct_str(100*kEff*sig_show(s))), 'Interpreter', 'latex');
    grid on; box on; set(gca, 'FontSize', fs);
end
exportgraphics(fig, fullfile(outFolder, 'slide_noise.png'), 'Resolution', 300);

%% ======================
% 3. Posterior variance vs noise (MCMC tuned vs Laplace) and correlation
% ======================
R = load(fullfile(P.results, 'compare_cov_noise_results_tuned.mat'));
sig = R.noise_levels;
sfine = logspace(log10(sig(1)), log10(sig(end)), 50);
Ew = R.exps(1);  Ep = R.exps(2);
col = lines(3);

fig = figure('Position', [100 100 1200 460]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
series = {Ew.std(1,:).^2, Ew.std(2,:).^2, Ep.std(2,:).^2};
lapl   = [Ew.std_unit(1), Ew.std_unit(2), Ep.std_unit(2)].^2;
lbl    = {'$\mathrm{Var}(a)$', '$\mathrm{Var}(w)$', '$\mathrm{Var}(p)$'};
% legend entries show marker (MCMC) + dashed line (Laplace) in the series colour
h = gobjects(1,3);
for k = 1:3
    loglog(sig, series{k}, 'o', 'MarkerSize', 8, 'LineWidth', 2, ...
        'Color', col(k,:), 'MarkerFaceColor', col(k,:)); hold on;
    loglog(sfine, lapl(k)*sfine.^2, '--', 'LineWidth', 1.5, 'Color', col(k,:));
    h(k) = loglog(nan, nan, 'o--', 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'Color', col(k,:), 'MarkerFaceColor', col(k,:));
end
% position is snapped to the grid: variance cannot drop below ~ (1/n)^2/12
hg = loglog(sig([1 end]), (1/n)^2/12*[1 1], ':', 'LineWidth', 2, 'Color', col(3,:));
noise_axis(gca, kEff);
xlabel('noise level $\sigma$ \ (effective relative noise)', 'Interpreter', 'latex');
ylabel('posterior variance', 'Interpreter', 'latex');
lg = legend([h hg], [lbl, {'grid limit for $p$'}], ...
    'Interpreter', 'latex', 'Location', 'southeast', 'NumColumns', 2);
title(lg, 'markers: MCMC, dashed: Laplace $\sigma^2 (G^{-1})_{ii}$', 'Interpreter', 'latex');
title('Variance $\propto \sigma^2$ (slope 2)', 'Interpreter', 'latex');
grid on; box on; set(gca, 'FontSize', fs);

nexttile;
semilogx(sig, Ew.corr, 'o', 'MarkerSize', 8, 'LineWidth', 2, ...
    'Color', col(1,:), 'MarkerFaceColor', col(1,:)); hold on;
semilogx(sig([1 end]), Ew.corr_lapl*[1 1], '--', 'LineWidth', 1.5, 'Color', col(1,:));
semilogx(sig, Ep.corr, 's', 'MarkerSize', 8, 'LineWidth', 2, ...
    'Color', col(2,:), 'MarkerFaceColor', col(2,:));
semilogx(sig([1 end]), Ep.corr_lapl*[1 1], '--', 'LineWidth', 1.5, 'Color', col(2,:));
h1 = semilogx(nan, nan, 'o--', 'MarkerSize', 8, 'LineWidth', 1.5, ...
    'Color', col(1,:), 'MarkerFaceColor', col(1,:));
h2 = semilogx(nan, nan, 's--', 'MarkerSize', 8, 'LineWidth', 1.5, ...
    'Color', col(2,:), 'MarkerFaceColor', col(2,:));
ylim([-1 1]);
noise_axis(gca, kEff);
xlabel('noise level $\sigma$ \ (effective relative noise)', 'Interpreter', 'latex');
ylabel('posterior correlation', 'Interpreter', 'latex');
lg = legend([h1 h2], {'$\mathrm{corr}(a,w)$', '$\mathrm{corr}(a,p)$'}, ...
    'Interpreter', 'latex', 'Location', 'southeast');
title(lg, 'markers: MCMC, dashed: Laplace', 'Interpreter', 'latex');
title('Correlation independent of $\sigma$', 'Interpreter', 'latex');
grid on; box on; set(gca, 'FontSize', fs);

exportgraphics(fig, fullfile(outFolder, 'slide_variance.png'), 'Resolution', 300);

% slopes of log Var vs log sigma (least squares, MCMC)
for k = 1:3
    pfit = polyfit(log(sig), log(series{k}), 1);
    fprintf('%s: fitted slope %.2f\n', lbl{k}, pfit(1));
end

%% ======================
% Local functions
% ======================
function [rho, x] = contrast(th, separation_mult, width_scale, C1_amp, C2_amp, n, h, lam, r, c0)
% rho(x) = a exp(-(x - x0(p))^4 / (2 w^2)), th = [a, w, p]
gen = MediumGenerator(th(3), separation_mult, width_scale, C1_amp, C2_amp);
gen.autoPlot = false;
[~, ~, x] = gen.medium3(n, h, lam, r, c0, th(1), th(2));
rho = gen.cont;
end
