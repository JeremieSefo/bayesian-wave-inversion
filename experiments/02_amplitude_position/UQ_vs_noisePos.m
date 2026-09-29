P = repo_paths();   % output folders (run after the MCMC script, which sets the workspace)
%% ============================================================
%  Noise vs Uncertainty study
%  Amplitude–Position inference
% ============================================================

% ============================================================
% Noise levels
% ============================================================

noise_levels = [1e-2, 5e-2, 1e-1, 2e-1, 4e-1, 8e-1];
nNoise = numel(noise_levels);

% ============================================================
% Storage
% ============================================================

amp_mean = zeros(nNoise,1);
amp_std  = zeros(nNoise,1);

pos_mean = zeros(nNoise,1);
pos_std  = zeros(nNoise,1);

amp_q025 = zeros(nNoise,1);
amp_q975 = zeros(nNoise,1);

pos_q025 = zeros(nNoise,1);
pos_q975 = zeros(nNoise,1);

corr_amp_pos = zeros(nNoise,1);

burnin_frac = 0.30;

% ============================================================
% Loop over noise levels
% ============================================================

for j = 1:nNoise

    sigma_noise = noise_levels(j);

    fprintf('Running noise level = %.3e\n', sigma_noise);

    % --------------------------------------------------------
    % Generate noisy observations
    % ---------------------------------------------------------

    dW = sqrt(dt) * randn(M, N-1);
    W  = [zeros(M,1), cumsum(dW,2)];

    X = X0 .* exp(TrueScattData_noisy.' + sigma_noise * W);

    % --------------------------------------------------------
    % Run MCMC
    % ---------------------------------------------------------

    theta = run_mcmc_wavePos( ...
        X, dt, sigma_noise, ...
        TrueOper, sigma_true, ...
        separation_mult, width_scale, ...
        C1_amp, C2_amp, ...
        n, hh, lam, r, c0);

    % --------------------------------------------------------
    % Remove burn-in
    % ---------------------------------------------------------

    burnin = round(burnin_frac * size(theta,2));

    amp_samp = theta(1, burnin+1:end);
    pos_samp = theta(2, burnin+1:end);

    % --------------------------------------------------------
    % Posterior summaries
    % ---------------------------------------------------------

    amp_mean(j) = mean(amp_samp);
    pos_mean(j) = mean(pos_samp);

    amp_std(j) = std(amp_samp);
    pos_std(j) = std(pos_samp);

    amp_q025(j) = quantile(amp_samp, 0.025);
    amp_q975(j) = quantile(amp_samp, 0.975);

    pos_q025(j) = quantile(pos_samp, 0.025);
    pos_q975(j) = quantile(pos_samp, 0.975);

    % --------------------------------------------------------
    % Posterior correlation
    % ---------------------------------------------------------

    R = corrcoef(amp_samp, pos_samp);
    corr_amp_pos(j) = R(1,2);

end

%% ============================================================
% 1. Posterior standard deviation vs noise
% ============================================================

fig = figure;

loglog(noise_levels, amp_std, 'o-', 'LineWidth', 1.5);
hold on;

loglog(noise_levels, pos_std, 's-', 'LineWidth', 1.5);

xlabel('Noise level $\sigma_{\mathrm{noise}}$', ...
    'Interpreter','latex');

ylabel('Posterior standard deviation', ...
    'Interpreter','latex');

legend({'Amplitude', 'Position'});

title('Posterior uncertainty scaling with noise', ...
    'Interpreter','latex');

grid on;

outFolder = P.report;
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

% saveTikz(fig,'uncertainty_vs_noise_position', outFolder)

%% ============================================================
% 2. Credible interval: amplitude
% ============================================================

fig = figure;

errorbar( ...
    noise_levels, ...
    amp_mean, ...
    amp_mean - amp_q025, ...
    amp_q975 - amp_mean, ...
    'o-', ...
    'LineWidth',1.5);

xlabel('$\sigma_{\mathrm{noise}}$', ...
    'Interpreter','latex');

ylabel('Amplitude $a$', ...
    'Interpreter','latex');

title('95\% credible interval for amplitude', ...
    'Interpreter','latex');

grid on;

% saveTikz(fig,'credible_interval_amp_vs_noise', outFolder)

%% ============================================================
% 3. Credible interval: position
% ============================================================

fig = figure;

errorbar( ...
    noise_levels, ...
    pos_mean, ...
    pos_mean - pos_q025, ...
    pos_q975 - pos_mean, ...
    'o-', ...
    'LineWidth',1.5);

xlabel('$\sigma_{\mathrm{noise}}$', ...
    'Interpreter','latex');

ylabel('Position $p$', ...
    'Interpreter','latex');

title('95\% credible interval for position', ...
    'Interpreter','latex');

grid on;

% saveTikz(fig,'credible_interval_pos_vs_noise', outFolder)

%% ============================================================
% 4. Posterior dependence structure
% ============================================================

fig = figure;

plot(noise_levels, corr_amp_pos, ...
    'o-', 'LineWidth',1.5);

xlabel('$\sigma_{\mathrm{noise}}$', ...
    'Interpreter','latex');

ylabel('Correlation$(a,p)$', ...
    'Interpreter','latex');

title('Posterior parameter dependence', ...
    'Interpreter','latex');

grid on;

% saveTikz(fig,'joint_dependence_vs_noise_position', outFolder)