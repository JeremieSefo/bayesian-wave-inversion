%% MCMC for wave speed estimation
clear; close all; clc;
rng(1);
%% load all variables from the MAT-file into workspace
P = repo_paths();
load(fullfile(P.data, 'waveData.mat'));
C1_true=C1;
% ======================
% 1. Generate synthetic data
% ======================
% Spatial grid
nI = 1500 * 2;
Nb  = 1;
n  = Nb*(nI-1) + 1;
Le = 2;
hh  = Le / (n-1);

% Frequencies
center_freq = 2*pi*0.02e3;
bandwidth   = 2*0.25*center_freq;
w = center_freq + bandwidth;
lambda = 2*pi / w;

pulse_index = 2;
c0 = 1;
f  = c0 / lambda;
lam = c0 / f;
r=2;

idx_frac1 = 7/50;
idx1 = idx_frac1 * n;
separation_mult = 1*3;
width_scale = 1;
C1_amp = 1. - .4;
C2_amp = 1. -0.;

amp_true = 0.75;
sigma_true = 0.5e-2;
sigNoise = 5e-3;

X0 = 1;
M=100;
T = 1.5;
N = numel(TrueScattData_noisy); % +9 4
dt = T/(N-1);
t = (0:N-1)*dt;

dW = sqrt(dt)*randn(M,N-1);
W = [zeros(M,1), cumsum(dW,2)];

X = X0 * exp(TrueScattData_noisy.' + sigNoise*W);
%% ======================
% 2. Plot a few trajectories
% ======================

nPlot = 10;                  % number of trajectories to display
idx = randperm(M, nPlot);   % randomly select paths

fig = figure;
hold on;

h = gobjects(nPlot,1);

for i = 1:nPlot
    h(i) = loglog(t, X(idx(i),:), 'LineWidth', 1);
end

h_det = loglog(t, X0*exp(TrueScattData_noisy.'), ...
    'k--', 'LineWidth', 2);

xlabel('t');
ylabel('X_t');
title('Sample trajectories of the measured scattered waves');
grid on;

legend([h; h_det], ...
    [arrayfun(@(i) sprintf('Path %d', i), 1:nPlot, 'UniformOutput', false), ...
     {'deterministic path'}], ...
    'Location', 'best');
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'trajectories', outFolder)
%% test 2 data sufficiency
negloglik = @(theta) -wave_loglik(theta(1), (theta(2)), X, dt,TrueOper, idx_frac1, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,hh,lam,r,c0,sigNoise);

theta0 = [.1; 1e-3];   % initial guess [mu, sigma]

theta_hat = fminsearch(negloglik, theta0);
%%
amp_hat = theta_hat(1);
sigma_hat = (theta_hat(2));

fprintf('MLE (numerical): amp = %.4f, sigma = %.4f\n', amp_hat, sigma_hat);

%% ======================
% 2. Log-likelihood function
% ======================
loglik = @(amp,sigma) wave_loglik(amp, sigma, X, dt,TrueOper, idx_frac1, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,hh,lam,r,c0, sigNoise);

% ======================
% 3. MCMC settings
% ======================
nIter = 4000;

theta = zeros(2,nIter);   % [amp; log(sigma)]
accept = 0;

% initial guess
theta(:,1) = [.1; 1e-3]; % [.75; 5e-3];

% proposal std
s_amp = 1e-2;
s_ls = 1e-3;

%% ======================
% 4. MCMC loop
% ======================
for i = 2:nIter
    
    amp_curr = theta(1,i-1);
    ls_curr = theta(2,i-1);
    sig_curr = (ls_curr); % exp
    
    logL_curr = loglik(amp_curr, sig_curr);
    
    % ---- proposal ----
    amp_prop = amp_curr + s_amp*randn;
    ls_prop = ls_curr + s_ls*randn;
    sig_prop = (ls_prop); % exp
    
    logL_prop = loglik(amp_prop, sig_prop);
    
    % ---- MH acceptance ----
    log_alpha = logL_prop - logL_curr;
    
    if log(rand) < log_alpha
        theta(:,i) = [amp_prop; ls_prop];
        accept = accept + 1;
    else
        theta(:,i) = theta(:,i-1);
    end
end

accept_rate = accept/nIter;
fprintf('Acceptance rate: %.2f%%\n', 100*accept_rate);

%% ======================
% MCMC error decay
% 5. Post-processing
% ======================
amp_samp = theta(1,:);
sig_samp = (theta(2,:)); % exp

burnin = round(0.3*nIter);
%MCMC estimates
amp_path=amp_samp(burnin:end);
sig_path=sig_samp(burnin:end);
amp_est = mean(amp_path);
sig_est = mean(sig_path);

fprintf('True amp = %.3f, Estimated amp = %.3f\n', amp_true, amp_est);
fprintf('True sigma = %.3f, Estimated sigma = %.3f\n', sigma_true, sig_est);
%%
amp_err = NaN(1,nIter);
sig_err = NaN(1,nIter);
burnin = 1; % round(0.3*nIter);
for k = burnin+100:nIter
    amp_est_k = mean(amp_samp(burnin:k));
    sig_est_k = mean(sig_samp(burnin:k));
    
    amp_err(k) = abs(amp_est_k - amp_true);
    sig_err(k) = abs(sig_est_k - sigma_true);
end
%%
idx = burnin+100:nIter;

fig=figure;

%fig = subplot(2,1,1);
loglog(idx, amp_err(idx), 'LineWidth', 1.5);
hold on;
ref = amp_err(idx(1)) * (idx/idx(1)).^(-0.5);
loglog(idx, ref, 'k--', 'LineWidth', 1);
legend('Error', 'O(N^{-1/2})');
xlabel('Number of MCMC samples');
ylabel('|a_{est} - a_{true}|');
title('Error decay of the amplitude');
grid on;
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'MCMC error deacy amp', outFolder)

fig=figure;
%fig=subplot(2,1,2);
loglog(idx, sig_err(idx), 'LineWidth', 1.5);
hold on;
ref = sig_err(idx(1)) * (idx/idx(1)).^(-0.5);
loglog(idx, ref, 'k--', 'LineWidth', 1);
legend('Error', 'O(N^{-1/2})');
xlabel('Number of MCMC samples');
ylabel('|w_{est} - w_{true}|');
title('Error decay of the width');
grid on;
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'MCMC error deacy sigma', outFolder)



%% ======================
% 6. MCMC path with update directions
% ======================

amp_path = theta(1,:);
sig_path = (theta(2,:)); % exp
burnin = 1;% round(0.1*nIter);
fig=figure;
hold on;

% plot chain path
plot(amp_path(burnin:end), sig_path(burnin:end), '*-', 'LineWidth', 1);

xlabel('amplitude a');
ylabel('width w');
title('MCMC update directions (chain evolution)');
grid on;
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'MCMC samples zoom', outFolder)
%% ======================
% 6. Plots
% ======================
fig=figure;
%subplot(2,1,1)
burnin =1;
plot(amp_path(burnin:end))
hold on
yline(amp_true,'r--')
title('amplitude chain')
grid on
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'amp chain', outFolder)

%subplot(2,1,2)
fig=figure;
plot(sig_path(burnin:end))
hold on
yline(sigma_true,'r--')
title('width chain')
grid on
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
saveTikz(fig,'sigma chain', outFolder)

fig=figure;
scatter(amp_path(burnin:end), sig_path(burnin:end), 25, '.')
xlabel('amplitude a'); ylabel('width w');
title('Posterior samples')
grid on
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'chain 2D', outFolder)

%ampWidth_uncerntainty;
%UQ_vs_noise;