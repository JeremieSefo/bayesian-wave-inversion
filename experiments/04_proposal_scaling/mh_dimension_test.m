%% How does random-walk MH scale with the number of parameters?
%  Target: exact Gaussian posterior of the Gaussian-basis coefficients
%  (sigma = 0.1, L2 prior with physical tau = 0.5, so the posterior is proper).
%  Samplers (optimal scaling 2.38^2/d):
%    diag : proposal N(0, (2.38^2/d) diag(Sigma))   (only marginal scales known)
%    full : proposal N(0, (2.38^2/d) Sigma)          (full posterior covariance known)
%  Reported: acceptance, min/median ESS over coordinates, ESS per 1000 iterations.
clear; clc;

P = repo_paths();
S = load(fullfile(P.data, 'waveData.mat'));
A = S.TrueOper;  D = S.TrueScattData_noisy;
[N, n] = size(A);  M = 100;  T = 1.5;  dt = T/(N-1);
x = linspace(0, 2, n).';  dx = x(2) - x(1);
J = diff(A, 1, 1);
sn = 0.1;  sn2 = sn^2*dt/M;  tau = 0.5;
rng(103);
dW = sqrt(dt)*randn(M, N-1);  W = [zeros(M,1), cumsum(dW,2)];
ybar = mean(diff(D.' + sn*W, 1, 2), 1).';

nIter = 40000;  burn = 5000;
fprintf('%4s %-5s | %6s %8s %8s %10s %7s | %s\n', 'd', 'prop', 'acc', 'minESS', 'medESS', 'ESS/1000it', 'time', 'cond(Sigma)');
for d = [2 20 80]
    nodes = linspace(0, 0.8, d);  h = nodes(2) - nodes(1);
    Phi = exp(-(x - nodes).^2 / (2*(0.6*h)^2));
    B = J*Phi;
    P = B.'*B/sn2 + (Phi.'*Phi)*dx/tau^2;  P = (P+P.')/2;
    Sigma = inv(P);  Sigma = (Sigma+Sigma.')/2;
    mu = P \ (B.'*ybar/sn2);
    logpost = @(be) -0.5*sum((ybar - B*be).^2)/sn2 - 0.5*(Phi*be).'*(Phi*be)*dx/tau^2;
    for prop = ["diag", "full"]
        if prop == "diag"
            L = sqrt((2.38^2/d) * diag(diag(Sigma)));
        else
            L = chol((2.38^2/d) * Sigma, 'lower');
        end
        rng(5);
        be = mu;  lp = logpost(be);  chain = zeros(d, nIter);  nAcc = 0;
        tic;
        for it = 1:nIter
            pr = be + L*randn(d,1);
            lq = logpost(pr);
            if log(rand) < lq - lp, be = pr; lp = lq; nAcc = nAcc + 1; end
            chain(:,it) = be;
        end
        tt = toc;
        smp = chain(:, burn+1:end);
        e = arrayfun(@(k) ess(smp(k,:)), 1:d);
        fprintf('%4d %-5s | %5.1f%% %8.0f %8.0f %10.2f %6.1fs | %.1e\n', d, prop, 100*nAcc/nIter, ...
            min(e), median(e), 1000*min(e)/size(smp,2), tt, cond(Sigma));
    end
end

function N_eff = ess(x)
x = x(:) - mean(x);  n = numel(x);  c0 = sum(x.^2);
if c0 == 0, N_eff = 1; return; end
tau = 1;
for k = 1:n-1
    r = sum(x(1:n-k).*x(1+k:n)) / c0;
    if r <= 0, break; end
    tau = tau + 2*r;
end
N_eff = n / tau;
end
