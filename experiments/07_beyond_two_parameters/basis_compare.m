%% Hat vs Gaussian basis for the Bayesian contrast inversion (flat prior)
%  For each basis and size Nb on [0, 0.8]:
%   - basis error      : best L2 fit of rho_true in the basis
%   - model error      : ||J (rho_true - Phi beta_best)|| relative to the noise
%                        norm of the averaged increments at sigma = 0.01
%   - posterior error  : ||c/c0 mean - true|| / ||true - 1|| on the bump region
%   - coverage         : fraction of the bump region inside the 95% band
%  Gaussian basis: phi_j(x) = exp(-(x - x_j)^2 / (2 (f h)^2)), h = node spacing.
clear; clc;

P = repo_paths();
S = load(fullfile(P.data, 'waveData.mat'));
A = S.TrueOper;  rho = S.Contrast;  D = S.TrueScattData_noisy;
[N, n] = size(A);  M = 100;  T = 1.5;  dt = T/(N-1);
x = linspace(0, 2, n).';
J = diff(A, 1, 1);
bump = x >= 0.1 & x <= 0.45;
ratio = @(r) 0.5 * (r + sqrt(r.^2 + 4));
ct = ratio(rho);

noise_levels = [1e-2, 5e-2, 1e-1, 2e-1, 4e-1, 8e-1];
jShow = [1 3 6];
NbList = [10 20 40 80 120];
f = 0.6;                                          % Gaussian width / spacing
noise_ref = sqrt(N-1) * 0.01 * sqrt(dt/M);        % noise norm at sigma = 0.01

% same averaged data as the other scripts
ybar = cell(1, numel(noise_levels));
for j = jShow
    rng(100 + j);
    dW = sqrt(dt)*randn(M, N-1);
    W  = [zeros(M,1), cumsum(dW,2)];
    ybar{j} = mean(diff(D.' + noise_levels(j)*W, 1, 2), 1).';
end

fprintf('%-6s %4s | %8s %9s %8s | %s\n', 'basis', 'Nb', 'basisErr', 'modelErr', 'cond', ...
    'posterior error / coverage at sigma = 0.01, 0.1, 0.8');
for basis = ["hat", "gauss"]
    for Nb = NbList
        nodes = linspace(0, 0.8, Nb);  h = nodes(2) - nodes(1);
        if basis == "hat"
            Phi = max(0, 1 - abs(x - nodes)/h);
        else
            Phi = exp(-(x - nodes).^2 / (2*(f*h)^2));
        end
        B = J*Phi;  BtB = B.'*B;
        bb = (Phi.'*Phi) \ (Phi.'*rho);
        berr = norm(Phi(bump,:)*bb - rho(bump)) / norm(rho(bump));
        merr = norm(J*(rho - Phi*bb)) / noise_ref;
        str = '';
        for j = jShow
            sn = noise_levels(j);
            mu = BtB \ (B.'*ybar{j});
            Sig = sn^2*dt/M * inv(BtB);  Sig = (Sig + Sig.')/2;
            [Vs, Ls] = eig(Sig);  Lc = Vs * sqrt(max(Ls, 0));      % robust square root
            rng(7);
            prof = ratio(Phi * (mu + Lc*randn(Nb, 2000)));
            lo = quantile(prof, 0.025, 2);  hi = quantile(prof, 0.975, 2);
            e  = norm(mean(prof(bump,:), 2) - ct(bump)) / norm(ct(bump) - 1);
            cv = mean(ct(bump) >= lo(bump) & ct(bump) <= hi(bump));
            str = [str sprintf('%7.1f%% /%4.0f%%  ', 100*e, 100*cv)]; %#ok<AGROW>
        end
        fprintf('%-6s %4d | %7.2f%% %8.0fx %8.1e | %s\n', basis, Nb, 100*berr, merr, cond(BtB), str);
    end
end
