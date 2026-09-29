%% Outlook figure: how many modes of the contrast can the data resolve?
%  J = diff(A) maps pixel values of rho to data increments. With the SVD
%  J = U S V', mode i of rho (coefficient c_i = v_i' rho) is measured with
%  posterior std  sigma / (s_i sqrt(M/dt)).  Mode i is resolvable when this
%  noise level is below |c_i| of the true medium (discrete Picard criterion).
clear; close all; clc;

P = repo_paths();
S = load(fullfile(P.data, 'waveData.mat'));
outFolder = P.slides;
fs = 16;

A = S.TrueOper;  rho = S.Contrast;
N = size(A, 1);  M = 100;  T = 1.5;  dt = T/(N-1);
kEff = rel_noise_factor(S.TrueScattData_noisy, T) / sqrt(M);

J = diff(A, 1, 1);
[~, Sg, V] = svd(J, 'econ');
s = diag(Sg);
c = abs(V' * rho);
i = (1:numel(s)).';

sig  = [0.01 0.1 0.8];
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];

fig = figure('Position', [100 100 1100 520]);
semilogy(i, c, '.', 'Color', [0.4 0.4 0.4], 'MarkerSize', 10); hold on;
cs = exp(movmean(log(c), 9));             % smoothed trend of |c_i|
semilogy(i, cs, 'k-', 'LineWidth', 2);
h = gobjects(1, numel(sig));
for k = 1:numel(sig)
    nl = sig(k) ./ (s * sqrt(M/dt));
    h(k) = semilogy(i, nl, '--', 'LineWidth', 2.2, 'Color', cols(k,:));
    nres = sum(nl < c);
    fprintf('sigma = %.2f (%s%%): %d resolvable modes\n', sig(k), pct_str(100*kEff*sig(k)), nres);
    lbl{k} = sprintf('noise, $\\sigma = %g$ (%s\\%%): %d modes', sig(k), pct_str(100*kEff*sig(k)), nres); %#ok<SAGROW>
end
xline(numel(s), 'k:', 'LineWidth', 1.5);
xlim([1 numel(s)+5]);
ylim([1e-6 1e3]);
xlabel('mode $i$ (singular vectors of $J = \Delta A$)', 'Interpreter', 'latex');
ylabel('magnitude', 'Interpreter', 'latex');
legend([findobj(gca, 'Type', 'line', 'LineStyle', 'none'); findobj(gca, 'Type', 'line', 'Color', [0 0 0], 'LineStyle', '-'); h(:)], ...
    [{'$|v_i^\top \rho_{\mathrm{true}}|$', 'trend'}, lbl], ...
    'Interpreter', 'latex', 'Location', 'northwest');
text(numel(s)-3, 2e-6, sprintf('%d data increments', numel(s)), 'HorizontalAlignment', 'right', ...
    'FontSize', fs-2, 'Interpreter', 'latex');
title('Modes of the contrast above the noise level', 'Interpreter', 'latex');
grid on; box on; set(gca, 'FontSize', fs);

exportgraphics(fig, fullfile(outFolder, 'slide_outlook.png'), 'Resolution', 300);
