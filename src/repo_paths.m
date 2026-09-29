function P = repo_paths()
%REPO_PATHS Folders of the repository, independent of where it is cloned.
%   P.root     repository root
%   P.data     input data (waveData.mat)
%   P.results  results saved by the experiments (.mat), reused by later scripts
%   P.figures  figures for the report (P.figures/report) and slides (P.slides)
%   P.slides   slide figures (PNG), with TikZ panels in P.slides/tikz
P.root    = fileparts(fileparts(mfilename('fullpath')));
P.data    = fullfile(P.root, 'data');
P.results = fullfile(P.root, 'results');
P.figures = fullfile(P.root, 'figures');
P.report  = fullfile(P.figures, 'report');
P.slides  = fullfile(P.figures, 'slides');
for f = {P.results, P.report, P.slides}
    if ~exist(f{1}, 'dir'), mkdir(f{1}); end
end
end
