%% SETUP  Add the repository folders to the MATLAB path (run once per session)
%  Run this from anywhere:  >> run('path/to/bayesian-wave-inversion/setup.m')
%  or open the file and press Run.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'data'));
exps = dir(fullfile(root, 'experiments', '0*'));
for k = 1:numel(exps)
    addpath(fullfile(exps(k).folder, exps(k).name));
end
P = repo_paths();   % also creates results/ and figures/ if missing
fprintf('Bayesian wave inversion: paths set (root = %s)\n', root);
if exist('matlab2tikz', 'file') ~= 2
    fprintf('  note: matlab2tikz not found; figures are saved as PNG only (TikZ export skipped)\n');
end
clear root exps k
