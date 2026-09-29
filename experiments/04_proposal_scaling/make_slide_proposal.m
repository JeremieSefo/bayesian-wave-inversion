%% Slide figure/table: original vs scaled proposal (amplitude-width)
%  original : fixed steps (1e-2, 1e-3), start (0.1, 1e-3), 30% burn-in
%             (compare_cov_noise_results.mat)
%  scaled   : production phase of the two-phase sampler, steps from G at the
%             best pilot state (slide_results.mat)
%  Both use the same data at each noise level and the same number of samples.
%  Reference: Laplace std and correlation from G at the true parameters
%  (compare_cov_noise_results_tuned.mat).
clear; close all; clc;

P = repo_paths();
O = load(fullfile(P.results, 'compare_cov_noise_results.mat'));
L = load(fullfile(P.results, 'compare_cov_noise_results_tuned.mat'));
S = load(fullfile(P.results, 'slide_results.mat'));
outFolder = P.slides;
fs = 15;
sig = O.noise_levels;
e = 1;                                        % amplitude-width
EO = O.exps(e);  EL = L.exps(e);
burnin = round(0.3 * size(EO.chains{1}, 2));

%% ======================
% Table: acceptance, ESS, std ratio, correlation
% ======================
fprintf('sigma | acc  ESS(a) std/Lap corr | acc  ESS(a) std/Lap corr   (Laplace corr %.2f)\n', EL.corr_lapl);
for j = 1:numel(sig)
    thO = EO.chains{j}(:, burnin+1:end);
    Ns  = size(thO, 2);
    thS = S.R(e,j).prod(:, 1:Ns);
    stdLap = sig(j) * EL.std_unit(1);
    rO = corrcoef(thO.');  rS = corrcoef(thS.');
    fprintf('%5.2f | %4.1f%% %5.0f %6.2f %+5.2f | %4.1f%% %5.0f %6.2f %+5.2f\n', sig(j), ...
        100*EO.acc(j), ess(thO(1,:)), std(thO(1,:))/stdLap, rO(1,2), ...
        100*S.R(e,j).accQ, ess(thS(1,:)), std(thS(1,:))/stdLap, rS(1,2));
end

%% ======================
% Figure: sigma = 0.01
% ======================
j  = 1;  sn = sig(j);
thO = EO.chains{j}(:, burnin+1:end);
Ns  = size(thO, 2);
thS = S.R(e,j).prod(:, 1:Ns);
it  = 1:Ns;
cO = [0 0.45 0.74];  cS = [0.85 0.33 0.10];
lblO = sprintf('original (acc. %.0f%%)', 100*EO.acc(j));
lblS = sprintf('scaled (acc. %.0f%%)',   100*S.R(e,j).accQ);

fig = figure('Position', [100 100 1300 760]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(it, thS(1,:), 'Color', cS, 'LineWidth', 0.7);
plot(it, thO(1,:), 'Color', cO, 'LineWidth', 2);
yline(0.75, 'k--', 'LineWidth', 1.2);
ylabel('amplitude $a$', 'Interpreter', 'latex');
xlabel('sample (after burn-in)', 'Interpreter', 'latex');
legend({lblS, lblO, 'true'}, 'Location', 'southeast');
title('Amplitude chain', 'Interpreter', 'latex');
ytickformat('%.3f'); grid on; box on; set(gca, 'FontSize', fs);

nexttile; hold on;
plot(it, 1e3*thS(2,:), 'Color', cS, 'LineWidth', 0.7);
plot(it, 1e3*thO(2,:), 'Color', cO, 'LineWidth', 2);
yline(5, 'k--', 'LineWidth', 1.2);
ylabel('width $w\ (\times 10^{-3})$', 'Interpreter', 'latex');
xlabel('sample (after burn-in)', 'Interpreter', 'latex');
title('Width chain', 'Interpreter', 'latex');
ytickformat('%.3f'); grid on; box on; set(gca, 'FontSize', fs);

% samples with the Laplace 95% ellipse (shape from G, centred at the posterior mean)
nexttile; hold on;
C = sn^2 * [EL.std_unit(1)^2, EL.cov_unit; EL.cov_unit, EL.std_unit(2)^2];
phi = linspace(0, 2*pi, 200);
m0 = mean(thS, 2);                       % Laplace posterior is centred at the estimate
ell = m0 + sqrtm(C) * sqrt(chi2inv(0.95, 2)) * [cos(phi); sin(phi)];
h1 = plot(thS(1,:), 1e3*thS(2,:), '.', 'Color', cS, 'MarkerSize', 5);
[uO, ~, icO] = unique(thO.', 'rows');
h2 = scatter(uO(:,1), 1e3*uO(:,2), 30 + 150*accumarray(icO,1)/Ns, cO, 'filled');
h3 = plot(ell(1,:), 1e3*ell(2,:), 'k-', 'LineWidth', 1.5);
xlabel('amplitude $a$', 'Interpreter', 'latex');
ylabel('width $w\ (\times 10^{-3})$', 'Interpreter', 'latex');
legend([h1 h2 h3], {'scaled', sprintf('original (%d states)', size(uO,1)), 'Laplace 95%'}, ...
    'Location', 'northwest');
title('Posterior samples', 'Interpreter', 'latex');
xtickformat('%.3f'); ytickformat('%.3f'); grid on; box on; set(gca, 'FontSize', fs);

nexttile; hold on;
maxLag = 200;
plot(0:maxLag, acf(thS(1,:), maxLag), 'Color', cS, 'LineWidth', 2);
plot(0:maxLag, acf(thO(1,:), maxLag), 'Color', cO, 'LineWidth', 2);
legend({sprintf('scaled (ESS = %.0f)', ess(thS(1,:))), ...
        sprintf('original (ESS = %.0f)', ess(thO(1,:)))}, 'Location', 'northeast');
xlabel('lag $j$', 'Interpreter', 'latex');
ylabel('autocorrelation $\rho_j$ of $a$', 'Interpreter', 'latex');
ylim([-0.2 1]);
title('Chain autocorrelation', 'Interpreter', 'latex');
grid on; box on; set(gca, 'FontSize', fs);

exportgraphics(fig, fullfile(outFolder, 'slide_proposal.png'), 'Resolution', 300);
fprintf('sigma = %.2f: original %d distinct states, %.1f%% inside Laplace 95%% ellipse\n', sn, size(uO,1), ...
    100*mean(sum((C\(thO - m0)).*(thO - m0), 1) <= chi2inv(0.95,2)));
fprintf('              scaled   %.1f%% inside Laplace 95%% ellipse\n', ...
    100*mean(sum((C\(thS - m0)).*(thS - m0), 1) <= chi2inv(0.95,2)));

%% ======================
% Local functions
% ======================
function r = acf(x, maxLag)
x = x(:) - mean(x);  n = numel(x);  c0 = sum(x.^2);
r = zeros(1, maxLag+1);
if c0 == 0, r(1) = 1; return; end
for k = 0:maxLag
    r(k+1) = sum(x(1:n-k) .* x(1+k:n)) / c0;
end
end

function N_eff = ess(x)
% n / tau, tau = 1 + 2 sum rho_k up to the first non-positive rho_k
n = numel(x);  r = acf(x, n-1);  tau = 1;
for k = 2:n
    if r(k) <= 0, break; end
    tau = tau + 2*r(k);
end
N_eff = n / tau;
end
