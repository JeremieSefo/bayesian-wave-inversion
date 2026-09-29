%% MAKE_DATA  Rebuild the true medium and the noise-free data from the operator
%  waveData.mat contains
%    TrueOper             270 x 3000  linear forward operator A: contrast -> scattered data
%    Contrast             3000 x 1    true contrast rho(x) on x in [0, 2]
%    C0, C1               3000 x 1    background and true wave speed (C1/C0 = c/c0)
%    TrueScattData_noisy  270 x 1     noise-free scattered data D = A*rho
%                                     (the name is historical: noise is added in the scripts)
%  Only TrueOper is primary data. This script regenerates Contrast, C0, C1 and
%  the data from the true parameters with MediumGenerator, checks them against
%  the stored file, and (optionally) writes a new waveData.mat.
%
%  To test another setting, change the parameters below (or the operator) and
%  set write = true. All experiments read data/waveData.mat.
clear; clc;
P = repo_paths();

% ---- true medium (as in the report and the slides) ----
amp   = 0.75;     % contrast amplitude a
width = 5e-3;     % width parameter w
pos   = 7/50;     % position p (fraction of the grid; x0 = 2p)
write = false;    % true: overwrite data/waveData.mat with the regenerated fields

% ---- grid and probing wave (fixed by the operator) ----
n  = 3000;  Le = 2;  h = Le/(n-1);
center_freq = 2*pi*0.02e3;
lam = 2*pi / (center_freq + 2*0.25*center_freq);

S = load(fullfile(P.data, 'waveData.mat'));
A = S.TrueOper;

gen = MediumGenerator(pos, 3, 1, 0.6, 1);
gen.autoPlot = false;
[C0, C1, x] = gen.medium3(n, h, lam, 2, 1, amp, width);
Contrast = gen.cont;
TrueScattData_noisy = A * Contrast;
TrueOper = A;

fprintf('regenerated vs stored (relative difference):\n');
fprintf('  Contrast %.1e | C0 %.1e | C1 %.1e | data %.1e\n', ...
    norm(Contrast - S.Contrast)/norm(S.Contrast), norm(C0 - S.C0)/norm(S.C0), ...
    norm(C1 - S.C1)/norm(S.C1), ...
    norm(TrueScattData_noisy - S.TrueScattData_noisy)/norm(S.TrueScattData_noisy));

if write
    save(fullfile(P.data, 'waveData.mat'), 'TrueOper', 'Contrast', 'C0', 'C1', ...
        'TrueScattData_noisy', '-v7');
    fprintf('written: %s\n', fullfile(P.data, 'waveData.mat'));
end

figure;
subplot(2,1,1); plot(x, C1./C0, 'LineWidth', 1.5); xlim([0 0.7]); grid on;
xlabel('x'); ylabel('c/c_0'); title('True medium');
subplot(2,1,2); plot(linspace(0, 1.5, size(A,1)), TrueScattData_noisy, 'LineWidth', 1.5); grid on;
xlabel('t'); ylabel('D^{sc}(t)'); title('Noise-free scattered data');
