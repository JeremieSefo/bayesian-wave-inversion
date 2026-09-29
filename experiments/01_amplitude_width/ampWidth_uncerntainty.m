P = repo_paths();   % output folders (run after the MCMC script, which sets the workspace)
burnin =  round(0.3*nIter);
amp_s = theta(1,burnin:end);
sig_s = theta(2,burnin:end);

Ns = numel(amp_s);
contrast_all = zeros(Ns, n);
for i = 1:Ns
    amp_i = amp_s(i);
    sig_i = sig_s(i);

    gennn = MediumGenerator(idx_frac1 , separation_mult, width_scale, ...
        C1_amp, C2_amp);
    gennn.autoPlot = false;

    [C00, C11, x] = gennn.medium3(n, hh, lam, r, c0, amp_i, sig_i);

    contrast_all(i,:) = C11 ./ C00;
end

contrast_mean = mean(contrast_all, 1);
contrast_low  = quantile(contrast_all, 0.025, 1);
contrast_high = quantile(contrast_all, 0.975, 1);
contrast_std = std(contrast_all, 0, 1);
%%
fig=figure; hold on;
nHalf = ceil(numel(x)/2);
fill([x(1:nHalf).' fliplr(x(1:nHalf).')], ...
     [contrast_low(1:nHalf) fliplr(contrast_high(1:nHalf))], ...
     [0.8 0.8 1], 'EdgeColor', 'none');

plot(x(1:nHalf), contrast_mean(1:nHalf), 'b', 'LineWidth', 1);
plot(x(1:nHalf)', C1_true(1:nHalf)./C0(1:nHalf), 'r', 'LineWidth', 1);
legend('Uncertainty region', 'Infered medium', 'True medium')
xlabel('x');
ylabel('C / C_0');
title('Posterior contrast field with uncertainty');
grid on;
outFolder = P.report; if ~exist(outFolder,'dir') mkdir(outFolder); end 
%saveTikz(fig,'contrast infered ', outFolder)