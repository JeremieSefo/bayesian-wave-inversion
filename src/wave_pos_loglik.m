function logL = wave_pos_loglik(amp, pos, X, dt, A, sigma, separation_mult, width_scale, ...
    C1_amp, C2_amp,n,h,lam,r,c0, sigNoise)
% X is M x (N+1)

Y = log(X);
dY = diff(Y,1,2);   % increments

% Generate medium (keep gennn local and avoid storing large fields)
gennn = MediumGenerator(pos, separation_mult, width_scale, ...
    C1_amp, C2_amp);
gennn.autoPlot = false;
gennn.medium3(n,h,lam,r,c0, amp, sigma);
rho=gennn.cont;
scatData = A*rho;
dscatData = diff(scatData,1);
var_inc  = sigNoise^2*dt;

logL = -0.5 * numel(dY) * log(2*pi*var_inc) ...
       - sum((dY - dscatData.').^2, 'all') / (2*var_inc);
end