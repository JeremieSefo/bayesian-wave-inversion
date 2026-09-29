P = repo_paths();   % output folders (run after the MCMC script, which sets the workspace)
%% ============================================================
%  Noise vs Uncertainty study for wave parameter inference
% ============================================================
% Noise levels to test
noise_levels = [1e-2, 5e-2, 1e-1, 2e-1, 4e-1, 8e-1];

nNoise = numel(noise_levels);

amp_mean = zeros(nNoise,1);
amp_std  = zeros(nNoise,1);
pos_mean = zeros(nNoise,1);
pos_std  = zeros(nNoise,1);

amp_q025 = zeros(nNoise,1);
amp_q975 = zeros(nNoise,1);

pos_q025 = zeros(nNoise,1);
pos_q975 = zeros(nNoise,1);

cov_amp_pos = zeros(nNoise,1);

burnin_frac = 0.3;

% ============================================================
% Loop over noise levels
% ============================================================
for j = 1:nNoise

    sigNoise = noise_levels(j);

    % --- generate noisy data ---
    dW = sqrt(dt)*randn(M,N-1);
    W = [zeros(M,1), cumsum(dW,2)];

    X = X0 * exp(TrueScattData_noisy.' + sigNoise*W);

    % --- run MCMC ---
    theta = run_mcmc_wave(X, dt, sigNoise,TrueOper, idx_frac1, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,hh,lam,r,c0);  % <-- wrap your code

    burnin = round(burnin_frac * size(theta,2));

    amp_samp = theta(1,burnin:end);
    pos_samp = theta(2,burnin:end);

    % ============================================================
    % Posterior summaries
    % ============================================================

    amp_mean(j) = mean(amp_samp);
    pos_mean(j) = mean(pos_samp);

    amp_std(j)  = std(amp_samp);
    pos_std(j)  = std(pos_samp);

    amp_q025(j) = quantile(amp_samp, 0.025);
    amp_q975(j) = quantile(amp_samp, 0.975);

    pos_q025(j) = quantile(pos_samp, 0.025);
    pos_q975(j) = quantile(pos_samp, 0.975);

    % joint uncertainty (covariance)
    C = cov(amp_samp, pos_samp);
    cov_amp_pos(j) = C(1,2);
end

%% ============================================================
% 1. Uncertainty scaling plots
% ============================================================

fig=figure;
loglog(noise_levels, amp_std, 'o-', 'LineWidth', 1.5); hold on;
loglog(noise_levels, pos_std, 's-', 'LineWidth', 1.5); hold on;
xlabel('Noise level \sigma_{noise}');
ylabel('Posterior standard deviation');
legend({'Amplitude', 'Width'});
title('Uncertainty scaling with noise');
grid on;
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'uncerntainty vs noise', outFolder)
%% ============================================================
% 2. Credible intervals
% ============================================================

fig=figure;
%subplot(1,2,1)
errorbar(noise_levels, amp_mean, ...
    amp_mean - amp_q025, amp_q975 - amp_mean, 'o-');
xlabel('\sigma_{noise}');
ylabel('amplitude a');
title('95% credible interval (amplitude)');
grid on;
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'credib interv amp vs noise', outFolder)

fig=figure;
%subplot(1,2,2)
errorbar(noise_levels, pos_mean, ...
    pos_mean - pos_q025, pos_q975 - pos_mean, 'o-');
xlabel('\sigma_{noise}');
ylabel('width w');
title('95% credible interval (Width)');
grid on;
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'credib interv width vs noise', outFolder)

%% ============================================================
% 3. Joint uncertainty structure
% ============================================================

fig=figure;
plot(noise_levels, cov_amp_pos, 'o-', 'LineWidth', 1.5);
xlabel('\sigma_{noise}');
ylabel('Covariance(amp, width)');
title('Joint posterior dependence');
grid on;
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'joint dependence vs noise', outFolder)