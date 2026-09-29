
%% SEE separate histograms  posterior densities
amp_path = theta(1,:);
sig_path = (theta(2,:));  % exp

burnin = 1;% round(0.1*nIter);

amp_samp = amp_path(burnin:end);
sig_samp = sig_path(burnin:end);

figure;

subplot(1,2,1)
histogram(amp_samp, 30, 'Normalization', 'pdf');
xline(amp_true, 'r--', 'LineWidth', 2);
xlabel('\amp');
ylabel('Density');
title('Posterior of \amp');
grid on;

subplot(1,2,2)
histogram(sig_samp, 30, 'Normalization', 'pdf');
xline(sigma_true, 'r--', 'LineWidth', 2);
xlabel('\sigma');
ylabel('Density');
title('Posterior of \sigma');
grid on;

%% SEE 2D histogram joint posterior density
amp_path = theta(1,:);
sig_path = (theta(2,:)); % exp

burnin = 1;%round(0.1*nIter);

amp_samp = amp_path(burnin:end);
sig_samp = sig_path(burnin:end);

figure;

histogram2(amp_samp, sig_samp, 40, ...
    'Normalization', 'pdf');

xlabel('\amp');
ylabel('\sigma');
title('Joint posterior distribution');
colorbar;
grid on;

%% SEE smooth Joint posterior density
figure;
burnin = 1500;% round(0.1*nIter);
%burnin = 1;
amp_samp = amp_path(burnin:end);
sig_samp = sig_path(burnin:end);
% grid for evaluation
amp_grid = linspace(min(amp_samp), max(amp_samp), 100);
sig_grid = linspace(min(sig_samp), max(sig_samp), 100);

[amp, SIG] = meshgrid(amp_grid, sig_grid);

% kernel density estimate
vals = [amp_samp(:), sig_samp(:)];
F = ksdensity(vals, [amp(:), SIG(:)]);
F = reshape(F, size(amp));

% contour plot
contourf(amp, SIG, F, 20, 'LineColor', 'none');
colormap(parula);   % or your preferred base map

c = colormap;
c(1,:) = [1 1 1];   % set lowest value (≈0) to white
colormap(c);
colorbar;

hold on;
scatter(amp_samp(1:5:end), sig_samp(1:5:end), 25, 'k', '.');

xlabel('\amp');
ylabel('\sigma');
title('Joint posterior');
grid on;