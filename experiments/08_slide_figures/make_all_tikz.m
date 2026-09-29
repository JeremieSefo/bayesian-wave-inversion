function make_all_tikz()
%MAKE_ALL_TIKZ Re-run all slide figure scripts and export their panels to TikZ
%  (figures/slides/tikz). A shim exportgraphics.m placed first on
%  the path replaces the PNG export, so the existing PNGs are not touched.
%  The scripts run in the base workspace, so their 'clear' does not affect
%  this function.
shim = fullfile(fileparts(mfilename('fullpath')), 'tikzshim');
P = repo_paths();
set(groot, 'defaultFigureVisible', 'off');
scripts = {'make_slide_figures', 'run_slide_results', 'make_slide_proposal', ...
           'make_slide_outlook', 'hat_basis_inversion', 'make_appendix_regularized', ...
           'make_slide_landscape', 'make_slide_ciwidth'};
addpath(shim, '-begin');
cleanup = onCleanup(@() rmpath(shim));
for i = 1:numel(scripts)
    fprintf('--- %s\n', scripts{i});
    evalin('base', sprintf('evalc(''%s'');', scripts{i}));
    close all;
    cd(P.root);
end
end
