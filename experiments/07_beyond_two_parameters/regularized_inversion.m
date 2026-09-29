%% Bayesian regularization on top of a basis expansion of the contrast
%  rho(x) = sum_j beta_j phi_j(x)  (hat or Gaussian functions on [0, 0.8])
%  Data   : ybar = B beta + e,  B = diff(A) Phi,  e ~ N(0, sn^2 I),  sn^2 = sigma^2 dt/M
%  Priors (Gaussian, zero mean), defined on the function rho:
%     'tikhonov' : ||rho||^2_{L2} / tau^2        -> P = Phi'Phi dx / tau^2
%     'gradient' : ||rho'||^2_{L2} / tau^2       -> P = dPhi'dPhi dx / tau^2  (+ tiny ridge)
%  Posterior: Sigma = (B'B/sn^2 + P)^{-1},  mu = Sigma B' ybar / sn^2
%  tau chosen by maximizing the evidence p(ybar | tau) (empirical Bayes).
%  Same noisy data as the other scripts (rng(100+j)).
clear; close all; clc;

P = repo_paths();
S = load(fullfile(P.data, 'waveData.mat'));
outFolder = P.slides;
fs = 13;

A = S.TrueOper;  rho = S.Contrast;  D = S.TrueScattData_noisy;
[N, n] = size(A);  M = 100;  T = 1.5;  dt = T/(N-1);
x = linspace(0, 2, n).';  dx = x(2) - x(1);
J = diff(A, 1, 1);
bump = x >= 0.1 & x <= 0.45;
ratio = @(r) 0.5 * (r + sqrt(r.^2 + 4));
ct = ratio(rho);
kEff = rel_noise_factor(D, T) / sqrt(M);

noise_levels = [1e-2, 5e-2, 1e-1, 2e-1, 4e-1, 8e-1];
jShow  = [1 3 6];
NbList = [40 80 120];
bases  = ["hat", "gauss"];
priors = ["flat", "tikhonov", "gradient", "tikfix", "gradfix"];
tauFix = struct('tikfix', 0.5, 'gradfix', 5);      % physical prior scales: ||rho|| ~ 0.3, ||rho''|| ~ 5
tauGrid = logspace(-4, 3, 57);
nDraw = 2000;

ybar = cell(1, numel(noise_levels));
for j = jShow
    rng(100 + j);
    dW = sqrt(dt)*randn(M, N-1);
    W  = [zeros(M,1), cumsum(dW,2)];
    ybar{j} = mean(diff(D.' + noise_levels(j)*W, 1, 2), 1).';
end

res = struct();
fprintf('%-6s %4s %-9s | %s\n', 'basis', 'Nb', 'prior', ...
    'err / coverage / tau at sigma = 0.01, 0.1, 0.8');
for basis = bases
    for Nb = NbList
        Phi = make_basis(basis, x, Nb);
        B = J*Phi;  BtB = B.'*B;
        dPhi = diff(Phi, 1, 1) / dx;
        K.tikhonov = (Phi.'*Phi) * dx;
        K.gradient = (dPhi.'*dPhi) * dx + 1e-8 * trace(dPhi.'*dPhi*dx)/Nb * eye(Nb);
        K.tikfix = K.tikhonov;  K.gradfix = K.gradient;
        for pr = priors
            str = '';
            for j = jShow
                sn2 = noise_levels(j)^2 * dt / M;
                y = ybar{j};
                if pr == "flat"
                    Pk = zeros(Nb);  tau = Inf;
                elseif isfield(tauFix, pr)
                    tau = tauFix.(pr);  Pk = K.(pr) / tau^2;
                else
                    % evidence: y ~ N(0, sn2 I + B Gamma B'),  Gamma = tau^2 K^{-1}
                    % B K^{-1} B' = G G' with G = B / chol(K); diagonalize once
                    Gk = B / chol((K.(pr) + K.(pr).')/2);
                    [Q, d] = eig(Gk*Gk.', 'vector');
                    d = max(d, 0);  z2 = (Q.'*y).^2;
                    logev = zeros(size(tauGrid));
                    for t = 1:numel(tauGrid)
                        ev = sn2 + tauGrid(t)^2 * d;          % eigenvalues of C
                        logev(t) = -0.5*sum(z2 ./ ev) - 0.5*sum(log(ev));
                    end
                    [~, it] = max(logev);  tau = tauGrid(it);
                    Pk = K.(pr) / tau^2;
                end
                Pp = BtB/sn2 + Pk;  Pp = (Pp + Pp.')/2;
                mu = Pp \ (B.'*y/sn2);
                [Vs, Ls] = eig(inv(Pp));  Lc = Vs * sqrt(max(Ls, 0));
                rng(7);
                prof = ratio(Phi * (mu + Lc*randn(Nb, nDraw)));
                lo = quantile(prof, 0.025, 2);  hi = quantile(prof, 0.975, 2);
                cm = mean(prof, 2);
                e  = norm(cm(bump) - ct(bump)) / norm(ct(bump) - 1);
                cv = mean(ct(bump) >= lo(bump) & ct(bump) <= hi(bump));
                key = sprintf('%s_%d_%s_%d', basis, Nb, pr, j);
                res.(key) = struct('cm', cm, 'lo', lo, 'hi', hi, 'err', e, 'cov', cv, 'tau', tau);
                str = [str sprintf('%6.1f%% /%4.0f%% /%8.1e   ', 100*e, 100*cv, tau)]; %#ok<AGROW>
            end
            fprintf('%-6s %4d %-9s | %s\n', basis, Nb, pr, str);
        end
    end
end
save(fullfile(P.results, 'regularized_results.mat'), 'res', 'noise_levels', 'NbList');

%% ======================
% Figure: Nb = 80, both bases, both priors, three noise levels
% ======================
Nb = 80;
rows = {"hat","tikhonov"; "hat","gradient"; "gauss","tikhonov"; "gauss","gradient"};
rowLbl = {'hat + $L^2$', 'hat + $\|\rho''\|^2$', 'Gauss + $L^2$', 'Gauss + $\|\rho''\|^2$'};
fig = figure('Position', [100 100 1400 1100]);
tiledlayout(size(rows,1), numel(jShow), 'TileSpacing', 'compact', 'Padding', 'compact');
xm = x <= 0.6;  xx = x(xm).';
for r = 1:size(rows,1)
    for q = 1:numel(jShow)
        j = jShow(q);
        R = res.(sprintf('%s_%d_%s_%d', rows{r,1}, Nb, rows{r,2}, j));
        nexttile; hold on;
        fill([xx fliplr(xx)], [R.lo(xm).' fliplr(R.hi(xm).')], [0.75 0.8 1], 'EdgeColor', 'none');
        plot(xx, R.cm(xm), 'b-', 'LineWidth', 1.6);
        plot(xx, ct(xm), 'r--', 'LineWidth', 1.4);
        xlim([0 0.6]);  ylim([0.7 1.8]);
        title(sprintf('%s, $\\sigma = %g$ (%s\\%%)', rowLbl{r}, noise_levels(j), ...
            pct_str(100*kEff*noise_levels(j))), 'Interpreter', 'latex');
        text(0.58, 1.72, sprintf('err %.0f\\%%, cov. %.0f\\%%', 100*R.err, 100*R.cov), ...
            'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'FontSize', fs-1);
        if q == 1, ylabel('$c/c_0$', 'Interpreter', 'latex'); end
        if r == size(rows,1), xlabel('$x$', 'Interpreter', 'latex'); end
        if r == 1 && q == 1
            legend({'95\% region', 'posterior mean', 'true'}, 'Interpreter', 'latex', ...
                'Location', 'west', 'FontSize', fs-3);
        end
        grid on; box on; set(gca, 'FontSize', fs);
    end
end
exportgraphics(fig, fullfile(outFolder, 'slide_regularized.png'), 'Resolution', 300);

%% ======================
% Local functions
% ======================
function Phi = make_basis(basis, x, Nb)
nodes = linspace(0, 0.8, Nb);  h = nodes(2) - nodes(1);
if basis == "hat"
    Phi = max(0, 1 - abs(x - nodes)/h);
else
    Phi = exp(-(x - nodes).^2 / (2*(0.6*h)^2));
end
end
