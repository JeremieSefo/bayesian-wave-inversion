function theta = run_mcmc_wavePos(X, dt, sigNoise, TrueOper, sigma, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,hh,lam,r,c0)

nIter = 8000;

theta = zeros(2,nIter);

theta(:,1) = [.1; 1/50];

s_amp = 1e-2;
s_pos = 1e-2;

for i = 2:nIter

    amp_curr = theta(1,i-1);
    pos_curr = theta(2,i-1);

    logL_curr = wave_pos_loglik(amp_curr, pos_curr, X, dt, ... 
        TrueOper, sigma, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,hh,lam,r,c0, sigNoise); 

    amp_prop = amp_curr + s_amp*randn;
    pos_prop = pos_curr + s_pos*randn;

    logL_prop = wave_pos_loglik(amp_prop, pos_prop, X, dt, ...
        TrueOper, sigma, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,hh,lam,r,c0, sigNoise);

    if log(rand) < (logL_prop - logL_curr)
        theta(:,i) = [amp_prop; pos_prop];
    else
        theta(:,i) = theta(:,i-1);
    end
end

end