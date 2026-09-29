%% Bayesian inversion of the contrast in a hat-function basis
%  rho(x) = sum_j beta_j phi_j(x), hat functions on [0, xmax] (illuminated region)
%  Data: mean log-increments ybar = B beta + noise, B = diff(A) * Phi,
%        noise ~ N(0, sigma^2 dt / M I)  (M trajectories averaged)
%  Flat prior => Gaussian posterior in closed form:
%        mu = (B'B)^{-1} B' ybar,   Sigma = sigma^2 dt/M (B'B)^{-1}
%  One random-walk MH chain checks that MCMC reproduces the exact posterior.
%  Same noisy data as compare_cov_noise.m (rng(100+j)).
clear; close all; clc;

P = repo_paths();
S = load(fullfile(P.data, 'waveData.mat'));
outFolder = P.slides;
fs = 14;

A = S.TrueOper;  rho_true = S.Contrast;  D = S.TrueScattData_noisy;
[N, n] = size(A);  M = 100;  T = 1.5;  dt = T/(N-1);
x = linspace(0, 2, n).';
J = diff(A, 1, 1);
kEff = rel_noise_factor(D, T) / sqrt(M);
ratio = @(r) 0.5 * (r + sqrt(r.^2 + 4));          % c/c0 from the contrast

noise_levels = [1e-2, 5e-2, 1e-1, 2e-1, 4e-1, 8e-1];
jShow  = [1 3 6];                                  % sigma = 0.01, 0.1, 0.8
NbList = [20 80 120];
xmax = 0.8;
nDraw = 2000;
bump = x >= 0.1 & x <= 0.45;

fig = figure('Position', [100 100 1400 900]);
tiledlayout(numel(NbList), numel(jShow), 'TileSpacing', 'compact', 'Padding', 'compact');
fprintf('%4s %6s | %8s %9s %9s %9s\n', 'Nb', 'sigma', 'cond B', 'err rho', 'err c/c0', 'coverage');
for b = 1:numel(NbList)
    Nb = NbList(b);
    nodes = linspace(0, xmax, Nb);
    hn = nodes(2) - nodes(1);
    Phi = max(0, 1 - abs(x - nodes) / hn);         % n x Nb hat functions
    B = J * Phi;
    BtB = B.' * B;
    beta_best = (Phi.'*Phi) \ (Phi.' * rho_true);  % best approximation in the basis
    for q = 1:numel(jShow)
        j = jShow(q);  sn = noise_levels(j);
        rng(100 + j);                               % same data as before
        dW = sqrt(dt)*randn(M, N-1);
        W  = [zeros(M,1), cumsum(dW,2)];
        X  = exp(D.' + sn*W);
        ybar = mean(diff(log(X), 1, 2), 1).';

        mu    = BtB \ (B.' * ybar);
        Sigma = sn^2 * dt / M * inv(BtB);
        Sigma = (Sigma + Sigma.') / 2;
        rng(7);
        betas = mu + chol(Sigma, 'lower') * randn(Nb, nDraw);
        prof  = ratio(Phi * betas);                % n x nDraw
        lo = quantile(prof, 0.025, 2);  hi = quantile(prof, 0.975, 2);
        cm = mean(prof, 2);
        ct = ratio(rho_true);

        err_rho = norm(Phi*mu - rho_true) / norm(rho_true);
        err_c   = norm(cm(bump) - ct(bump)) / norm(ct(bump) - 1);
        cover   = mean(ct(bump) >= lo(bump) & ct(bump) <= hi(bump));
        fprintf('%4d %6.2f | %8.1e %8.1f%% %8.1f%% %8.1f%%\n', Nb, sn, cond(BtB), ...
            100*err_rho, 100*err_c, 100*cover);

        nexttile; hold on;
        xm = x <= 0.6;  xx = x(xm).';
        fill([xx fliplr(xx)], [lo(xm).' fliplr(hi(xm).')], [0.75 0.8 1], 'EdgeColor', 'none');
        plot(xx, cm(xm), 'b-', 'LineWidth', 1.6);
        plot(xx, ct(xm), 'r--', 'LineWidth', 1.4);
        plot(xx, ratio(Phi(xm,:)*beta_best), 'k:', 'LineWidth', 1.2);
        xlim([0 0.6]);  ylim([0.5 2]);
        title(sprintf('$N_b = %d$, $\\sigma = %g$ (%s\\%%)', Nb, sn, pct_str(100*kEff*sn)), ...
            'Interpreter', 'latex');
        if q == 1, ylabel('$c/c_0$', 'Interpreter', 'latex'); end
        if b == numel(NbList), xlabel('$x$', 'Interpreter', 'latex'); end
        if b == 1 && q == 1
            legend({'95\% region', 'posterior mean', 'true', 'best fit in basis'}, ...
                'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fs-3);
        end
        grid on; box on; set(gca, 'FontSize', fs);
    end
end
exportgraphics(fig, fullfile(outFolder, 'slide_hat.png'), 'Resolution', 300);

%% ======================
% MCMC check: Nb = 20, sigma = 0.1
% ======================
Nb = 20;  j = 3;  sn = noise_levels(j);
nodes = linspace(0, xmax, Nb);  hn = nodes(2) - nodes(1);
Phi = max(0, 1 - abs(x - nodes) / hn);
B = J * Phi;  BtB = B.'*B;
rng(100 + j);
dW = sqrt(dt)*randn(M, N-1);  W = [zeros(M,1), cumsum(dW,2)];
ybar = mean(diff(log(exp(D.' + sn*W)), 1, 2), 1).';
mu = BtB \ (B.'*ybar);  Sigma = sn^2*dt/M * inv(BtB);  Sigma = (Sigma+Sigma.')/2;

% random-walk MH with full proposal covariance (2.38^2/d) Sigma, start at 0
var_inc = sn^2 * dt / M;
loglik = @(be) -sum((ybar - B*be).^2) / (2*var_inc);
Lp = chol((2.38^2/Nb) * Sigma, 'lower');
nIter = 60000;  be = zeros(Nb,1);  ll = loglik(be);  chain = zeros(Nb, nIter);  nAcc = 0;
rng(11);
for it = 1:nIter
    prop = be + Lp*randn(Nb,1);
    lp = loglik(prop);
    if log(rand) < lp - ll, be = prop; ll = lp; nAcc = nAcc + 1; end
    chain(:,it) = be;
end
smp = chain(:, 20001:end);
fprintf('\nMCMC check (Nb = %d, sigma = %.2f): acceptance %.1f%%\n', Nb, sn, 100*nAcc/nIter);
fprintf('  max |mean_MCMC - mu| / std  = %.2f\n', max(abs(mean(smp,2) - mu) ./ sqrt(diag(Sigma))));
fprintf('  std_MCMC / std_exact: min %.2f, max %.2f\n', min(std(smp,0,2)./sqrt(diag(Sigma))), ...
    max(std(smp,0,2)./sqrt(diag(Sigma))));
