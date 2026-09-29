%% Numerical results for the slides: inversion (chains, estimates) and UQ
%  Two-phase Metropolis-Hastings for amplitude-width and amplitude-position:
%   1. pilot phase from a far initial guess with fixed steps (as in the report)
%   2. production phase with steps s_i = 1.7 sigma sqrt((G^{-1})_ii), where G
%      is evaluated at the best pilot state (no use of the true parameters)
%  Posterior summaries and uncertainty regions use the production samples.
%  Same noisy data as compare_cov_noise.m (rng(100+j)).
clear; close all; clc;

P = repo_paths();
load(fullfile(P.data, 'waveData.mat'));
outFolder = P.slides;
if ~exist(outFolder,'dir') mkdir(outFolder); end
fs = 15;

%% ======================
% Setup (as in WAVE_MCMC / WAVEpos_MCMC)
% ======================
nI = 1500 * 2;  n = nI;  Le = 2;  hh = Le / (n-1);
center_freq = 2*pi*0.02e3;
w   = center_freq + 2*0.25*center_freq;
lam = 2*pi / w;
c0 = 1;  r = 2;
idx_frac1 = 7/50;
separation_mult = 3;  width_scale = 1;
C1_amp = 0.6;  C2_amp = 1;
amp_true = 0.75;  width_true = 5e-3;  pos_true = idx_frac1;

X0 = 1;  M = 100;  T = 1.5;
N  = numel(TrueScattData_noisy);
dt = T/(N-1);
kEff = rel_noise_factor(TrueScattData_noisy, T) / sqrt(M);   % effective relative noise per unit sigma
noise_levels = [1e-2, 5e-2, 1e-1, 2e-1, 4e-1, 8e-1];
nNoise = numel(noise_levels);

nPilot = 3000;       % pilot (burn-in) iterations
nProd  = 5000;       % production iterations

% rho for theta = [a; w] or [a; p]
rhoW = @(th) medium_rho(th(1), th(2), idx_frac1, separation_mult, width_scale, C1_amp, C2_amp, n, hh, lam, r, c0);
rhoP = @(th) medium_rho(th(1), width_true, th(2), separation_mult, width_scale, C1_amp, C2_amp, n, hh, lam, r, c0);

exps(1) = struct('name', 'width', 'sym', 'w', 'rho', rhoW, 'th_true', [amp_true; width_true], ...
    'th0', [.1; 1e-3], 's0', [1e-2; 1e-3], 'hmin', [0; 0]);
exps(2) = struct('name', 'position', 'sym', 'p', 'rho', rhoP, 'th_true', [amp_true; pos_true], ...
    'th0', [.1; 1/50], 's0', [1e-2; 1e-2], 'hmin', [0; 1/n]);

%% ======================
% Two-phase MCMC for every noise level
% ======================
for e = 1:2
    fwd = @(th) TrueOper * exps(e).rho(th);
    for j = 1:nNoise
        sn = noise_levels(j);
        rng(100 + j);                               % same data as compare_cov_noise.m
        dW = sqrt(dt)*randn(M, N-1);
        W  = [zeros(M,1), cumsum(dW,2)];
        X  = X0 * exp(TrueScattData_noisy.' + sn*W);
        dY = diff(log(X), 1, 2);

        rng(200 + 10*e + j);                        % sampler randomness
        % phase 1: pilot from a far initial guess, fixed steps
        [thP, logLP, accP] = run_mh(fwd, dY, sn, dt, exps(e).th0, exps(e).s0, nPilot);
        [~, ib] = max(logLP);
        th_best = thP(:, ib);

        % phase 2: steps from the Laplace approximation at the best pilot state
        J = jacobian_increments(fwd, th_best, exps(e).hmin);
        Cunit = inv(J.'*J) / M * dt;                % Sigma_post / sigma^2
        s = max(1.7 * sn * sqrt(diag(Cunit)), exps(e).hmin);
        [thQ, ~, accQ] = run_mh(fwd, dY, sn, dt, th_best, s, nProd);

        R(e,j).pilot = thP;   R(e,j).prod = thQ;
        R(e,j).accP  = accP;  R(e,j).accQ = accQ;  R(e,j).s = s;
        R(e,j).mean  = mean(thQ, 2);
        R(e,j).std   = std(thQ, 0, 2);
        R(e,j).q025  = quantile(thQ, 0.025, 2);
        R(e,j).q975  = quantile(thQ, 0.975, 2);
        fprintf('%-8s sigma %.0e: acc pilot %4.1f%%, prod %4.1f%% | a = %.4f [%.4f, %.4f], %s = %.5f [%.5f, %.5f]\n', ...
            exps(e).name, sn, 100*accP, 100*accQ, R(e,j).mean(1), R(e,j).q025(1), R(e,j).q975(1), ...
            exps(e).sym, R(e,j).mean(2), R(e,j).q025(2), R(e,j).q975(2));
    end
end
save(fullfile(P.results, 'slide_results.mat'), 'R', 'exps', 'noise_levels', 'nPilot', 'nProd');

%% ======================
% Figure 1: Markov chains at the lowest noise level
% ======================
j = 1;
fig = figure('Position', [100 100 1200 640]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ylab = {'amplitude $a$', 'width $w$'; 'amplitude $a$', 'position $p$'};
for e = 1:2
    th = [R(e,j).pilot, R(e,j).prod];
    it = 1:size(th, 2);
    for k = 1:2
        nexttile; hold on;
        plot(it(1:nPilot), th(k,1:nPilot), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
        plot(it(nPilot+1:end), th(k,nPilot+1:end), 'Color', [0 0.45 0.74], 'LineWidth', 0.8);
        yline(exps(e).th_true(k), 'r--', 'LineWidth', 1.5);
        xline(nPilot, 'k:', 'LineWidth', 1.5);
        ylabel(ylab{e,k}, 'Interpreter', 'latex');
        if e == 2, xlabel('MCMC iteration', 'Interpreter', 'latex'); end
        if e == 1 && k == 1
            legend({'pilot (burn-in)', 'production', 'true value'}, 'Location', 'southeast');
        end
        grid on; box on; set(gca, 'FontSize', fs);
    end
end
exportgraphics(fig, fullfile(outFolder, 'slide_chains.png'), 'Resolution', 300);

%% ======================
% Figure 2: uncertainty region of c/c0 for both experiments
% ======================
jShow = [1 3 5 6];                 % sigma = 0.01, 0.1, 0.4, 0.8
[~, x] = rhoW(exps(1).th_true);
xm = x <= 0.5;
ratio = @(rho) 0.5 * (rho + sqrt(rho.^2 + 4));      % c/c0 from the contrast
thin = 5;

fig = figure('Position', [100 100 1400 620]);
tiledlayout(2, numel(jShow), 'TileSpacing', 'compact', 'Padding', 'compact');
for e = 1:2
    rtrue = ratio(exps(e).rho(exps(e).th_true));
    for q = 1:numel(jShow)
        j = jShow(q);
        smp = R(e,j).prod(:, 1:thin:end);
        prof = zeros(size(smp,2), nnz(xm));
        for i = 1:size(smp,2)
            ri = ratio(exps(e).rho(smp(:,i)));
            prof(i,:) = ri(xm);
        end
        lo = quantile(prof, 0.025, 1);
        hi = quantile(prof, 0.975, 1);
        mu = mean(prof, 1);
        xx = x(xm).';
        nexttile; hold on;
        hf = fill([xx fliplr(xx)], [lo fliplr(hi)], [0.75 0.8 1], 'EdgeColor', 'none');
        hm = plot(xx, mu, 'b-', 'LineWidth', 1.8);
        ht = plot(xx, rtrue(xm), 'r--', 'LineWidth', 1.5);
        xlim([0.1 0.45]); ylim([0.98 1.75]);
        title(sprintf('$\\sigma = %g$ \\ (%s\\%%)', noise_levels(j), ...
            pct_str(100*kEff*noise_levels(j))), 'Interpreter', 'latex');
        if q == 1
            ylabel({sprintf('$(a,%s)$', exps(e).sym), '$c/c_0$'}, 'Interpreter', 'latex');
        end
        if e == 2, xlabel('$x$', 'Interpreter', 'latex'); end
        if e == 1 && q == 1
            legend([hf hm ht], {'95\% region', 'inferred', 'true'}, ...
                'Interpreter', 'latex', 'Location', 'northeast');
        end
        grid on; box on; set(gca, 'FontSize', fs);
    end
end
exportgraphics(fig, fullfile(outFolder, 'slide_uq_region.png'), 'Resolution', 300);

%% ======================
% Figure 3: posterior means and 95% credible intervals vs noise
% ======================
fig = figure('Position', [100 100 1400 430]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
panels = {1, 1, 'amplitude $a$ \ $(a,w)$'; 1, 2, 'width $w$'; 2, 2, 'position $p$'};
for k = 1:3
    e = panels{k,1};  i = panels{k,2};
    m  = arrayfun(@(jj) R(e,jj).mean(i), 1:nNoise);
    lo = arrayfun(@(jj) R(e,jj).q025(i), 1:nNoise);
    hi = arrayfun(@(jj) R(e,jj).q975(i), 1:nNoise);
    nexttile; hold on;
    errorbar(noise_levels, m, m - lo, hi - m, 'o-', 'LineWidth', 1.8, ...
        'MarkerSize', 7, 'MarkerFaceColor', [0 0.45 0.74], 'CapSize', 8);
    yline(exps(e).th_true(i), 'r--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log');
    xlim([0.007 1.1]);
    noise_axis(gca, kEff);
    xlabel('noise level $\sigma$ \ (effective relative noise)', 'Interpreter', 'latex');
    title(panels{k,3}, 'Interpreter', 'latex');
    if k == 1
        legend({'mean $\pm$ 95\% CI', 'true'}, 'Interpreter', 'latex', 'Location', 'northwest');
    end
    ytickformat('%.4g');
    grid on; box on; set(gca, 'FontSize', fs);
end
exportgraphics(fig, fullfile(outFolder, 'slide_ci.png'), 'Resolution', 300);

%% ======================
% Local functions
% ======================
function [rho, x] = medium_rho(amp, width, pos, separation_mult, width_scale, C1_amp, C2_amp, n, h, lam, r, c0)
gen = MediumGenerator(pos, separation_mult, width_scale, C1_amp, C2_amp);
gen.autoPlot = false;
[~, ~, x] = gen.medium3(n, h, lam, r, c0, amp, width);
rho = gen.cont;
end

function J = jacobian_increments(fwd, th, hmin)
J = zeros(numel(fwd(th))-1, numel(th));
for k = 1:numel(th)
    dth = zeros(size(th));
    dth(k) = max(1e-4 * abs(th(k)), hmin(k));
    J(:,k) = diff(fwd(th + dth) - fwd(th - dth)) / (2*dth(k));
end
end

function [theta, logL, acc] = run_mh(fwd, dY, sigNoise, dt, th0, s, nIter)
% Random-walk Metropolis-Hastings, flat prior
var_inc = sigNoise^2 * dt;
loglik  = @(th) -sum((dY - diff(fwd(th)).').^2, 'all') / (2*var_inc);
theta = zeros(2, nIter);  logL = zeros(1, nIter);
theta(:,1) = th0;  logL(1) = loglik(th0);
nAcc = 0;
for i = 2:nIter
    th_prop = theta(:,i-1) + s .* randn(2,1);
    lp = loglik(th_prop);
    if log(rand) < lp - logL(i-1)
        theta(:,i) = th_prop;  logL(i) = lp;  nAcc = nAcc + 1;
    else
        theta(:,i) = theta(:,i-1);  logL(i) = logL(i-1);
    end
end
acc = nAcc / (nIter-1);
end
