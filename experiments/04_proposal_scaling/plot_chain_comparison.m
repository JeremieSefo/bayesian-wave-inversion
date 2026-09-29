%% Report figures: original vs tuned Metropolis-Hastings chains
%  Uses the results saved by compare_cov_noise.m (run it with
%  useTuned = false and useTuned = true first).
clear; close all; clc;

P = repo_paths();
orig  = load(fullfile(P.results, 'compare_cov_noise_results.mat'));
tuned = load(fullfile(P.results, 'compare_cov_noise_results_tuned.mat'));
noise_levels = tuned.noise_levels;
burnin_frac = 0.3;

outFolder = fullfile(P.report, 'chain_compare');
if ~exist(outFolder,'dir') mkdir(outFolder); end

%% ======================
% 1. Acceptance rate and effective sample size, all noise levels
% ======================
names = {'width', 'position'};
fprintf('%-9s %7s | %8s %8s %8s | %8s %8s %8s\n', '', 'sigma', ...
    'acc orig', 'ESS a', 'ESS lam', 'acc tun', 'ESS a', 'ESS lam');
for e = 1:2
    for j = 1:numel(noise_levels)
        [ea_o, el_o] = chain_ess(orig.exps(e).chains{j}, burnin_frac);
        [ea_t, el_t] = chain_ess(tuned.exps(e).chains{j}, burnin_frac);
        fprintf('%-9s %7.0e | %7.1f%% %8.0f %8.0f | %7.1f%% %8.0f %8.0f\n', names{e}, ...
            noise_levels(j), 100*orig.exps(e).acc(j), ea_o, el_o, ...
            100*tuned.exps(e).acc(j), ea_t, el_t);
    end
end

%% ======================
% 2. Chains for amplitude-width at the lowest noise level
% ======================
j  = 1;                                   % sigma_noise = 1e-2
sn = noise_levels(j);
Eo = orig.exps(1);  Et = tuned.exps(1);
thO = Eo.chains{j};  thT = Et.chains{j};
nIter = size(thO, 2);
b   = round(burnin_frac * nIter);
it  = b+1:nIter;
wsc = 1e3;                                % plot width in units of 1e-3

[essAo, essWo] = chain_ess(thO, burnin_frac);
[essAt, essWt] = chain_ess(thT, burnin_frac);
lblO = sprintf('original (acc. %.0f%%)', 100*Eo.acc(j));
lblT = sprintf('tuned (acc. %.0f%%)',    100*Et.acc(j));

% ---- trace of the amplitude ----
fig = figure;
plot(it, thO(1,it), 'LineWidth', 1.2); hold on;
plot(it, thT(1,it), 'LineWidth', 0.8);
xlabel('MCMC iteration', 'Interpreter', 'latex');
ylabel('amplitude $a$', 'Interpreter', 'latex');
legend({lblO, lblT}, 'Location', 'southeast');
title('Amplitude chain', 'Interpreter', 'latex');
grid on;
saveTikz(fig, 'trace_amp', outFolder)

% ---- trace of the width ----
fig = figure;
plot(it, wsc*thO(2,it), 'LineWidth', 1.2); hold on;
plot(it, wsc*thT(2,it), 'LineWidth', 0.8);
xlabel('MCMC iteration', 'Interpreter', 'latex');
ylabel('width $w \; (\times 10^{-3})$', 'Interpreter', 'latex');
legend({lblO, lblT}, 'Location', 'southeast');
title('Width chain', 'Interpreter', 'latex');
grid on;
saveTikz(fig, 'trace_width', outFolder)

% ---- posterior samples with Laplace 95% ellipse ----
Cunit = [Et.std_unit(1)^2, Et.cov_unit; Et.cov_unit, Et.std_unit(2)^2];
C  = sn^2 * Cunit;
m  = mean(thT(:,it), 2);
phi = linspace(0, 2*pi, 200);
ell = m + sqrtm(C) * sqrt(chi2inv(0.95, 2)) * [cos(phi); sin(phi)];

[uO, ~, icO] = unique(thO(:,it).', 'rows');   % distinct states of the original chain
nO = accumarray(icO, 1);

fig = figure;
plot(thT(1,it), wsc*thT(2,it), '.', 'MarkerSize', 4, 'Color', [0.85 0.33 0.10]); hold on;
scatter(uO(:,1), wsc*uO(:,2), 20 + 80*nO/max(nO), [0 0.45 0.74], 'filled');
plot(ell(1,:), wsc*ell(2,:), 'k-', 'LineWidth', 1.2);
xlabel('amplitude $a$', 'Interpreter', 'latex');
ylabel('width $w \; (\times 10^{-3})$', 'Interpreter', 'latex');
legend({'tuned', sprintf('original (%d distinct states)', size(uO,1)), ...
    'Laplace 95% region'}, 'Location', 'northwest');
title('Posterior samples', 'Interpreter', 'latex');
grid on;
saveTikz(fig, 'samples_2d', outFolder)

% ---- autocorrelation of the amplitude chain ----
maxLag = 200;
fig = figure;
plot(0:maxLag, acf(thO(1,it), maxLag), 'LineWidth', 1.5); hold on;
plot(0:maxLag, acf(thT(1,it), maxLag), 'LineWidth', 1.5);
xlabel('lag', 'Interpreter', 'latex');
ylabel('autocorrelation of $a$', 'Interpreter', 'latex');
legend({sprintf('original (ESS = %.0f)', essAo), sprintf('tuned (ESS = %.0f)', essAt)}, ...
    'Location', 'northeast');
title('Chain autocorrelation', 'Interpreter', 'latex');
ylim([-0.2 1]);
grid on;
saveTikz(fig, 'acf_amp', outFolder)

fprintf('\nsigma = %.0e, amplitude-width, %d samples after burn-in:\n', sn, numel(it));
fprintf('  original: acc %.1f%%, ESS(a) = %.0f, ESS(w) = %.0f, distinct states = %d\n', ...
    100*Eo.acc(j), essAo, essWo, size(uO,1));
fprintf('  tuned   : acc %.1f%%, ESS(a) = %.0f, ESS(w) = %.0f\n', 100*Et.acc(j), essAt, essWt);
fprintf('  proposal steps: original [%.1e %.1e], tuned [%.1e %.1e]\n', Eo.s, ...
    1.7*sn*Et.std_unit);

%% ======================
% Local functions
% ======================
function r = acf(x, maxLag)
% Sample autocorrelation of x for lags 0..maxLag
x = x(:) - mean(x);
n = numel(x);
c0 = sum(x.^2);
r = zeros(1, maxLag+1);
if c0 == 0, r(1) = 1; return; end
for k = 0:maxLag
    r(k+1) = sum(x(1:n-k) .* x(1+k:n)) / c0;
end
end

function [essA, essL] = chain_ess(theta, burnin_frac)
% Effective sample size of both parameters after burn-in
b = round(burnin_frac * size(theta,2));
essA = ess(theta(1, b+1:end));
essL = ess(theta(2, b+1:end));
end

function N_eff = ess(x)
% ESS = n / tau, tau = 1 + 2 sum_k rho_k, truncated at the first
% non-positive autocorrelation (initial positive sequence)
n = numel(x);
r = acf(x, n-1);
tau = 1;
for k = 2:n
    if r(k) <= 0, break; end
    tau = tau + 2*r(k);
end
N_eff = n / tau;
end
