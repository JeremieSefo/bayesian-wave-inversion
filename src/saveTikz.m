function [] = saveTikz(fig, fileName, outdir, xlabels)
%SAVE TIKZ and PNG of the specified figure. Optionally set custom x tick labels for TikZ.
% Usage: saveTikz(fig, fileName, outdir, xlabels)
% xlabels should be a cell array of strings or empty to leave unchanged.

if nargin < 3 || isempty(outdir), outdir = 'texPlots'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
if nargin < 4, xlabels = []; end

outfile = fullfile(outdir, strcat(fileName, '.tex'));
pngfile = fullfile(outdir, strcat(fileName, '.png'));

% Validate figure handle
if ~ishandle(fig) || ~strcmp(get(fig, 'Type'), 'figure')
    error('First argument must be a valid figure handle.');
end

% Save PNG of specified figure
print(fig, pngfile, '-dpng', '-r300');

% Apply custom x tick labels if provided
if ~isempty(xlabels)
    ax = findall(fig, 'type', 'axes');
    % If multiple axes, apply to each visible axes that has XTick
    for k = 1:numel(ax)
        t = get(ax(k), 'XTick');
        if ~isempty(t)
            % Ensure xlabels length matches number of ticks or is scalar/cell for formatting
            if iscell(xlabels) && numel(xlabels) == numel(t)
                set(ax(k), 'XTickLabel', xlabels);
            elseif ischar(xlabels) || (iscell(xlabels) && isscalar(xlabels))
                % Single format string or single label: replicate
                set(ax(k), 'XTickLabel', repmat({xlabels}, size(t)));
            else
                warning('Provided xlabels do not match number of XTicks; skipping XTickLabel assignment for an axis.');
            end
        end
    end
end

% TikZ export needs matlab2tikz (optional dependency); the PNG is always saved
if exist('matlab2tikz', 'file') ~= 2
    fprintf('saveTikz: matlab2tikz not found, saved %s only\n', pngfile);
    return
end
figure(fig);        % make target figure current
cleanfigure();
matlab2tikz(outfile, ...
    'width','\figurewidth', ...
    'height','\figureheight', ...
    'interpretTickLabels', false);
end
function [] = saveTikz0(fig, fileName, outdir)
%SAVE TIKZ and PNG of the specified figure.
if nargin < 3 || isempty(outdir), outdir = 'texPlots'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
outfile = fullfile(outdir, strcat(fileName, '.tex'));
pngfile = fullfile(outdir, strcat(fileName, '.png'));
% Save PNG of specified figure
if ~ishandle(fig)
    error('First argument must be a valid figure handle.');
end
print(gcf, pngfile, '-dpng', '-r300');
% TikZ export needs matlab2tikz (optional dependency); the PNG is always saved
if exist('matlab2tikz', 'file') ~= 2
    fprintf('saveTikz: matlab2tikz not found, saved %s only\n', pngfile);
    return
end
figure(fig);        % make target figure current
cleanfigure();
matlab2tikz(outfile, ...
    'width','\figurewidth', ...
    'height','\figureheight', ...
    'interpretTickLabels', false);
end