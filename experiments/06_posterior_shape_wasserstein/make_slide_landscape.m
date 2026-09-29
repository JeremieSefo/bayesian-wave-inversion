%% Appendix figure: misfit landscape in the position p (a, w fixed at true values)
%  Left : L2 misfit of the increments (= -log L), noise-free and sigma = 0.1,
%         normalized and on a log scale to show the local minima
%  Right: W2^2 between positive densities exp(c s) of the mean log-signal and
%         the model signal, for noise-free data and sigma = 0.01, 0.1
%  Ybar(t) = mean over the M = 100 trajectories of log X_t (same data, rng(103)).
clear; close all; clc;

P = repo_paths();
S = load(fullfile(P.data, 'waveData.mat'));
outFolder = P.slides;
fs = 15;

A = S.TrueOper;  D = S.TrueScattData_noisy;
[N, n] = size(A);  M = 100;  T = 1.5;  dt = T/(N-1);  hh = 2/(n-1);
t = (0:N-1).' * dt;
cf = 2*pi*20;  lam = 2*pi/(cf + 0.5*cf);
p0 = 7/50;
kEff = rel_noise_factor(D, T) / sqrt(M);

pg = linspace(0.06, 0.30, 481);
d = zeros(N, numel(pg));
for k = 1:numel(pg)
    d(:,k) = A * medium_rho(0.75, 5e-3, pg(k), n, hh, lam);
end

% mean log-signal for each noise level (sigma = 0: noise-free)
sig = [0 0.01 0.1];
Yb = zeros(N, numel(sig));
for q = 1:numel(sig)
    rng(103);
    dW = sqrt(dt)*randn(M, N-1);  W = [zeros(M,1), cumsum(dW,2)];
    Yb(:,q) = mean(D.' + sig(q)*W, 1).';
end

c = 3 / max(abs(D));
pos = @(s) exp(c*s);                                  % positive density transform
L2 = @(y) sum((diff(y) - diff(d)).^2, 1);
W2 = @(y) arrayfun(@(k) w2sq(pos(y), pos(d(:,k)), t), 1:numel(pg));
nrm = @(v) (v - min(v)) / (max(v) - min(v));
isMin = @(v) find(v(2:end-1) < v(1:end-2) & v(2:end-1) < v(3:end)) + 1;

cols = [0 0 0; 0 0.45 0.74; 0.85 0.33 0.10];
fig = figure('Position', [100 100 1300 500]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
lab = {};  h = gobjects(0);
for q = [1 3]
    v = nrm(L2(Yb(:,q))) + 1e-6;
    h(end+1) = semilogy(pg, v, 'LineWidth', 2, 'Color', cols(q,:)); %#ok<SAGROW>
    m = isMin(v);
    semilogy(pg(m), v(m), 'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', cols(q,:));
    lab{end+1} = sprintf('%s: %d local minima', noise_lbl(sig(q), kEff), numel(m)); %#ok<SAGROW>
end
set(gca, 'YScale', 'log');
xline(p0, 'r--', 'LineWidth', 1.5);
xlabel('position $p$', 'Interpreter', 'latex');
ylabel('normalized $-\log L$ (log scale)', 'Interpreter', 'latex');
title('$L^2$ misfit (likelihood)', 'Interpreter', 'latex');
legend(h, lab, 'Interpreter', 'latex', 'Location', 'southeast');
grid on; box on; set(gca, 'FontSize', fs);

nexttile; hold on;
lab = {};  h = gobjects(0);
for q = 1:3
    v = nrm(W2(Yb(:,q)));
    h(end+1) = plot(pg, v, 'LineWidth', 2, 'Color', cols(q,:)); %#ok<SAGROW>
    m = isMin(v);
    plot(pg(m), v(m), 'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', cols(q,:));
    [~, i] = min(v);
    lab{end+1} = sprintf('%s: %d min., at $p = %.3f$', noise_lbl(sig(q), kEff), numel(m), pg(i)); %#ok<SAGROW>
end
xline(p0, 'r--', 'LineWidth', 1.5);
xlabel('position $p$', 'Interpreter', 'latex');
ylabel('normalized $W_2^2$', 'Interpreter', 'latex');
title('Wasserstein misfit', 'Interpreter', 'latex');
legend(h, lab, 'Interpreter', 'latex', 'Location', 'north');
grid on; box on; set(gca, 'FontSize', fs);

exportgraphics(fig, fullfile(outFolder, 'slide_landscape.png'), 'Resolution', 300);

%% ======================
% Local functions
% ======================
function s = noise_lbl(sn, kEff)
if sn == 0, s = 'noise-free';
else, s = sprintf('$\\sigma = %g$ (%s\\%%)', sn, pct_str(100*kEff*sn)); end
end

function v = w2sq(f, g, t)
% 1-D Wasserstein-2 distance squared between normalized densities f, g on t
f = f(:) / sum(f);  g = g(:) / sum(g);
F = cumsum(f);  G = cumsum(g);
s = linspace(0, 1, 2001).';  s = s(2:end-1);
[Fu, iF] = unique(F);  [Gu, iG] = unique(G);
v = mean((interp1(Fu, t(iF), s, 'linear', 'extrap') - interp1(Gu, t(iG), s, 'linear', 'extrap')).^2);
end

function rho = medium_rho(a, w, p, n, h, lam)
gen = MediumGenerator(p, 3, 1, 0.6, 1);
gen.autoPlot = false;
gen.medium3(n, h, lam, 2, 1, a, w);
rho = gen.cont;
end
