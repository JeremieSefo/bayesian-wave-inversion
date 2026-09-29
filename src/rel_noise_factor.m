function k = rel_noise_factor(D, T)
%REL_NOISE_FACTOR Relative noise level per unit sigma for one trajectory.
%   k = E[ RMS(W) ] / RMS(D), so that ||sigma W|| / ||D|| ~ k * sigma, where
%   D is the noise-free scattered data on [0, T] and W a standard Brownian
%   motion on the same grid. With M independent trajectories the effective
%   relative noise is k * sigma / sqrt(M).
N  = numel(D);
dt = T / (N-1);
rs = RandStream('mt19937ar', 'Seed', 0);        % leaves the global rng untouched
nMC = 2000;
W  = [zeros(nMC,1), cumsum(sqrt(dt) * randn(rs, nMC, N-1), 2)];
k  = mean(sqrt(mean(W.^2, 2))) / sqrt(mean(D(:).^2));
end
