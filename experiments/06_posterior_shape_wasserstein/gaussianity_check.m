%% Is the posterior Gaussian away from the true parameters?
%  1. Profile of log L along w and along p (a fixed at its true value),
%     compared with the Laplace quadratic from G at the truth.
%  2. Skewness / excess kurtosis of the production MCMC samples (slide_results.mat).
clear; clc;

P = repo_paths();
S = load(fullfile(P.data, 'waveData.mat'));
A = S.TrueOper;  D = S.TrueScattData_noisy;
[N, n] = size(A);  M = 100;  T = 1.5;  dt = T/(N-1);  hh = 2/(n-1);
cf = 2*pi*20;  lam = 2*pi/(cf + 0.5*cf);
rhoFun = @(a, w, p) medium_rho(a, w, p, n, hh, lam);
a0 = 0.75;  w0 = 5e-3;  p0 = 7/50;

sn = 0.1;
rng(103);
dW = sqrt(dt)*randn(M, N-1);  W = [zeros(M,1), cumsum(dW,2)];
ybar = mean(diff(D.' + sn*W, 1, 2), 1).';
sn2 = sn^2*dt/M;
loglik = @(rho) -sum((ybar - diff(A*rho)).^2) / (2*sn2);

% ---- profile along w ----
wg = linspace(1e-3, 1.2e-2, 400);
lw = arrayfun(@(w) loglik(rhoFun(a0, w, p0)), wg);
% ---- profile along p ----
pg = linspace(0.03, 0.30, 800);
lp = arrayfun(@(p) loglik(rhoFun(a0, w0, p)), pg);

% Laplace curvatures (1-D, a fixed) by finite differences at the peak
[~, iw] = max(lw);  [~, ip] = max(lp);
cw = -(lw(iw+1) - 2*lw(iw) + lw(iw-1)) / (wg(2)-wg(1))^2;
cp_ = -(lp(ip+1) - 2*lp(ip) + lp(ip-1)) / (pg(2)-pg(1))^2;
sw = 1/sqrt(cw);  sp = 1/sqrt(cp_);
fprintf('sigma = %.2f: 1-D posterior std (a fixed): w %.2e, p %.2e\n', sn, sw, sp);

% where does the quadratic stop describing the profile? (log-L error > 1)
qw = lw(iw) - 0.5*cw*(wg - wg(iw)).^2;   qp = lp(ip) - 0.5*cp_*(pg - pg(ip)).^2;
bw = abs(lw - qw) > 1;  bp = abs(lp - qp) > 1;
fprintf('w: quadratic accurate (error < 1 in log L) for |w - w_hat| < %.1f std\n', ...
    min(abs(wg(bw) - wg(iw))) / sw);
fprintf('p: quadratic accurate (error < 1 in log L) for |p - p_hat| < %.1f std\n', ...
    min(abs(pg(bp) - pg(ip))) / sp);

% local maxima of the profiles
isMax = @(v) find(v(2:end-1) > v(1:end-2) & v(2:end-1) > v(3:end)) + 1;
mw = isMax(lw);  mp = isMax(lp);
fprintf('w profile: %d local maxima on [%.0e, %.0e]\n', numel(mw), wg(1), wg(end));
fprintf('p profile: %d local maxima on [%.2f, %.2f]; the 3 highest secondary ones (p, drop in log L):\n', ...
    numel(mp), pg(1), pg(end));
[vals, order] = sort(lp(mp), 'descend');
for k = 2:min(4, numel(order))
    fprintf('   p = %.3f  (%.1f std from peak), log L lower by %.0f\n', pg(mp(order(k))), ...
        abs(pg(mp(order(k))) - pg(ip))/sp, lp(ip) - vals(k));
end

% ---- moments of the MCMC samples ----
R = load(fullfile(P.results, 'slide_results.mat'));
names = {'a (a,w)', 'w', 'a (a,p)', 'p'};
fprintf('\nMCMC production samples: skewness / excess kurtosis (Gaussian: 0 / 0)\n');
for j = [1 3 6]
    str = sprintf('sigma = %.2f: ', R.noise_levels(j));
    k = 0;
    for e = 1:2
        th = R.R(e,j).prod;
        for i = 1:2
            k = k + 1;
            str = [str sprintf('%s %+.2f/%+.2f   ', names{k}, skewness(th(i,:)), kurtosis(th(i,:)) - 3)]; %#ok<AGROW>
        end
    end
    fprintf('%s\n', str);
end

function rho = medium_rho(a, w, p, n, h, lam)
gen = MediumGenerator(p, 3, 1, 0.6, 1);
gen.autoPlot = false;
gen.medium3(n, h, lam, 2, 1, a, w);
rho = gen.cont;
end
