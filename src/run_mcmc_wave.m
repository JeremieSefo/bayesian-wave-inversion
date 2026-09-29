function theta = run_mcmc_wave(X, dt, sigNoise, TrueOper, idx_frac1, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,hh,lam,r,c0)

nIter = 4000;

theta = zeros(2,nIter);

theta(:,1) = [.1; 1e-3];

s_amp = 1e-2;
s_pos = 1e-3;

for i = 2:nIter

    amp_curr = theta(1,i-1);
    pos_curr = theta(2,i-1);

    logL_curr = wave_loglik(amp_curr, pos_curr, X, dt, ... 
        TrueOper, idx_frac1, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,hh,lam,r,c0, sigNoise); 

    amp_prop = amp_curr + s_amp*randn;
    pos_prop = pos_curr + s_pos*randn;

    logL_prop = wave_loglik(amp_prop, pos_prop, X, dt, ...
        TrueOper, idx_frac1, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,hh,lam,r,c0, sigNoise);

    if log(rand) < (logL_prop - logL_curr)
        theta(:,i) = [amp_prop; pos_prop];
    else
        theta(:,i) = theta(:,i-1);
    end
end

end