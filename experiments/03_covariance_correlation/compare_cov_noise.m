%% Covariance / correlation vs noise: amplitude-width and amplitude-position
%  Combines WAVE_MCMC + UQ_vs_noise and WAVEpos_MCMC + UQ_vs_noisePos in
%  one script, using the same noisy data for both experiments at each noise
%  level. Reports posterior covariance AND correlation, and compares them
%  with the small-noise (Laplace/Fisher) prediction
%       Sigma_post ~ inv(F),  F = J'J * M / (sigma^2 dt),
%  where J is the Jacobian of the data increments w.r.t. (a, lambda).
clearvars -except useTuned; close all; clc;
rng(1);

%% ======================
% 1. Setup (as in WAVE_MCMC / WAVEpos_MCMC)
% ======================
P = repo_paths();
load(fullfile(P.data, 'waveData.mat'));

nI = 1500 * 2;
Nb = 1;
n  = Nb*(nI-1) + 1;
Le = 2;
hh = Le / (n-1);

center_freq = 2*pi*0.02e3;
bandwidth   = 2*0.25*center_freq;
w = center_freq + bandwidth;
lambda = 2*pi / w;
c0  = 1;
f   = c0 / lambda;
lam = c0 / f;
r   = 2;

idx_frac1 = 7/50;
separation_mult = 1*3;
width_scale = 1;
C1_amp = 1. - .4;
C2_amp = 1. - 0.;

amp_true   = 0.75;
width_true = 0.5e-2;
pos_true   = idx_frac1;

X0 = 1;
M  = 100;
T  = 1.5;
N  = numel(TrueScattData_noisy);
dt = T/(N-1);

noise_levels = [1e-2, 5e-2, 1e-1, 2e-1, 4e-1, 8e-1];
nNoise = numel(noise_levels);
burnin_frac = 0.3;

% Sampler settings.
%   useTuned = false : original settings (reproduces Figs. 6 and 9)
%   useTuned = true  : start at the truth, proposal scaled to the Laplace std
if ~exist('useTuned', 'var'), useTuned = false; end

% Forward maps theta -> scattered data (N x 1)
fwdW = @(th) forward(th(1), th(2), idx_frac1, TrueOper, separation_mult, width_scale, ...
    C1_amp, C2_amp, n, hh, lam, r, c0);          % theta = [a; w]
fwdP = @(th) forward(th(1), width_true, th(2), TrueOper, separation_mult, width_scale, ...
    C1_amp, C2_amp, n, hh, lam, r, c0);          % theta = [a; p]

exps(1) = struct('name', 'width',    'fwd', fwdW, 'th_true', [amp_true; width_true], ...
    'th0', [.1; 1e-3], 's', [1e-2; 1e-3], 'nIter', 4000, 'hmin', [0; 0]);
exps(2) = struct('name', 'position', 'fwd', fwdP, 'th_true', [amp_true; pos_true], ...
    'th0', [.1; 1/50], 's', [1e-2; 1e-2], 'nIter', 8000, 'hmin', [0; 1/n]);
% NB: MediumGenerator snaps the position to the grid (idx1 = floor(n*p)+1),
% so the likelihood is piecewise constant in p with steps of 1/n. 'hmin' is
% the smallest meaningful step per parameter.

%% ======================
% 2. Laplace (Fisher) prediction at the true parameters
% ======================
for e = 1:2
    J = jacobian_increments(exps(e).fwd, exps(e).th_true, exps(e).hmin);
    exps(e).JtJ = J.'*J;                        % F = JtJ * M / (sigma^2 dt)
    Cunit = inv(exps(e).JtJ) / M * dt;          % covariance / sigma^2
    exps(e).cov_unit  = Cunit(1,2);
    exps(e).std_unit  = sqrt(diag(Cunit));
    exps(e).corr_lapl = Cunit(1,2) / sqrt(Cunit(1,1)*Cunit(2,2));
    fprintf('Laplace (%s): corr(a,%s) = %+.3f, Cov/sigma^2 = %+.3e\n', ...
        exps(e).name, exps(e).name, exps(e).corr_lapl, exps(e).cov_unit);
end

%% ======================
% 3. MCMC for each noise level, same data for both experiments
% ======================
for e = 1:2
    exps(e).mean = zeros(2,nNoise);  exps(e).std = zeros(2,nNoise);
    exps(e).q025 = zeros(2,nNoise);  exps(e).q975 = zeros(2,nNoise);
    exps(e).cov  = zeros(1,nNoise);  exps(e).corr = zeros(1,nNoise);
    exps(e).acc  = zeros(1,nNoise);  exps(e).chains = cell(1,nNoise);
end

for j = 1:nNoise
    sigNoise = noise_levels(j);

    rng(100 + j);                                % same data for both experiments
    dW = sqrt(dt)*randn(M,N-1);
    W  = [zeros(M,1), cumsum(dW,2)];
    X  = X0 * exp(TrueScattData_noisy.' + sigNoise*W);
    dY = diff(log(X), 1, 2);

    for e = 1:2
        if useTuned
            th0 = exps(e).th_true;
            s   = max(1.7 * sigNoise * exps(e).std_unit, exps(e).hmin);
        else
            th0 = exps(e).th0;
            s   = exps(e).s;
        end

        fprintf('noise %.0e, %-8s ... ', sigNoise, exps(e).name);
        tic;
        [theta, acc] = run_mh(exps(e).fwd, dY, sigNoise, dt, th0, s, exps(e).nIter);
        fprintf('acc %4.1f%%, %.0f s\n', 100*acc, toc);

        burnin = round(burnin_frac * exps(e).nIter);
        smp = theta(:, burnin+1:end);

        C = cov(smp.');
        exps(e).mean(:,j) = mean(smp, 2);
        exps(e).std(:,j)  = std(smp, 0, 2);
        exps(e).q025(:,j) = quantile(smp, 0.025, 2);
        exps(e).q975(:,j) = quantile(smp, 0.975, 2);
        exps(e).cov(j)    = C(1,2);
        exps(e).corr(j)   = C(1,2) / sqrt(C(1,1)*C(2,2));
        exps(e).acc(j)    = acc;
        exps(e).chains{j} = theta;
    end
end

resfile = 'compare_cov_noise_results.mat';
if useTuned, resfile = 'compare_cov_noise_results_tuned.mat'; end
save(fullfile(P.results, resfile), 'exps', 'noise_levels', 'useTuned');

%% ======================
% 4. Summary table
% ======================
for e = 1:2
    fprintf('\n=== amplitude-%s ===\n', exps(e).name);
    fprintf('%8s %12s %12s %9s %9s %7s\n', 'sigma', 'Cov MCMC', 'Cov Laplace', ...
        'corr', 'corr Lap', 'acc');
    for j = 1:nNoise
        fprintf('%8.0e %+12.3e %+12.3e %+9.3f %+9.3f %6.1f%%\n', noise_levels(j), ...
            exps(e).cov(j), noise_levels(j)^2*exps(e).cov_unit, ...
            exps(e).corr(j), exps(e).corr_lapl, 100*exps(e).acc(j));
    end
end

%% ======================
% 5. Plots
% ======================
for e = 1:2
    fig = figure('Name', ['amplitude-' exps(e).name]);

    subplot(1,2,1);
    plot(noise_levels, exps(e).cov, 'o-', 'LineWidth', 1.5); hold on;
    plot(noise_levels, noise_levels.^2 * exps(e).cov_unit, 'k--', 'LineWidth', 1);
    yline(0, ':');
    xlabel('\sigma_{noise}'); ylabel(['Cov(a, ' exps(e).name ')']);
    legend('MCMC', 'Laplace \sigma^2 (F^{-1})_{12}', 'Location', 'best');
    title('Posterior covariance'); grid on;

    subplot(1,2,2);
    plot(noise_levels, exps(e).corr, 'o-', 'LineWidth', 1.5); hold on;
    yline(exps(e).corr_lapl, 'k--', 'LineWidth', 1);
    yline(0, ':');
    ylim([-1 1]);
    xlabel('\sigma_{noise}'); ylabel(['corr(a, ' exps(e).name ')']);
    legend('MCMC', 'Laplace', 'Location', 'best');
    title('Posterior correlation'); grid on;
end

% Posterior samples at each noise level (after burn-in)
for e = 1:2
    figure('Name', ['samples amplitude-' exps(e).name]);
    for j = 1:nNoise
        subplot(2,3,j);
        th = exps(e).chains{j};
        burnin = round(burnin_frac * exps(e).nIter);
        plot(th(1,burnin+1:end), th(2,burnin+1:end), '.');
        xlabel('a'); ylabel(exps(e).name);
        title(sprintf('\\sigma = %.0e, corr = %+.2f', noise_levels(j), exps(e).corr(j)));
        grid on;
    end
end

%% ======================
% Local functions
% ======================
function d = forward(amp, width, pos, A, separation_mult, width_scale, ...
    C1_amp, C2_amp, n, h, lam, r, c0)
% Scattered data A*rho for rho = amp*exp(-(x-pos)^4/(2 width^2))
gen = MediumGenerator(pos, separation_mult, width_scale, C1_amp, C2_amp);
gen.autoPlot = false;
gen.medium3(n, h, lam, r, c0, amp, width);
d = A * gen.cont;
end

function J = jacobian_increments(fwd, th, hmin)
% Central-difference Jacobian of the data increments diff(fwd(th));
% step at least hmin (one grid cell for the position)
J = zeros(numel(fwd(th))-1, numel(th));
for k = 1:numel(th)
    dth = zeros(size(th));
    dth(k) = max(1e-4 * abs(th(k)), hmin(k));
    J(:,k) = diff(fwd(th + dth) - fwd(th - dth)) / (2*dth(k));
end
end

function [theta, acc] = run_mh(fwd, dY, sigNoise, dt, th0, s, nIter)
% Random-walk Metropolis-Hastings with flat prior (same likelihood as
% wave_loglik / wave_pos_loglik, constant terms dropped)
var_inc = sigNoise^2 * dt;
loglik  = @(th) -sum((dY - diff(fwd(th)).').^2, 'all') / (2*var_inc);

theta = zeros(2, nIter);
theta(:,1) = th0;
logL_curr = loglik(th0);
nAcc = 0;
for i = 2:nIter
    th_prop = theta(:,i-1) + s .* randn(2,1);
    logL_prop = loglik(th_prop);
    if log(rand) < logL_prop - logL_curr
        theta(:,i) = th_prop;
        logL_curr = logL_prop;
        nAcc = nAcc + 1;
    else
        theta(:,i) = theta(:,i-1);
    end
end
acc = nAcc / (nIter-1);
end
