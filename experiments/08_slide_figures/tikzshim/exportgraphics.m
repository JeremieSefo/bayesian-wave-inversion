function exportgraphics(fig, file, varargin)
% Shim used only by make_all_tikz.m: instead of writing the PNG, export the
% figure panels to TikZ in <png folder>/tikz/<name>_<k>.tex
[folder, name] = fileparts(file);
export_tikz_panels(fig, name, fullfile(folder, 'tikz'));
end
