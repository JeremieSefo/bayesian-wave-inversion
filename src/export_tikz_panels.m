function export_tikz_panels(fig, base, outFolder)
%EXPORT_TIKZ_PANELS Export every panel of a (tiled) figure to its own TikZ file.
%   Panels are written as <base>_<k>.tex (k = tile index, row-major) in
%   outFolder, with width \figurewidth and height \figureheight so that the
%   LaTeX side sets the size. matlab2tikz does not know tiledlayout or
%   xline/yline, so each axes (with its legend) is copied into a fresh figure
%   and constant lines are converted to ordinary lines first.
if exist('matlab2tikz', 'file') ~= 2
    fprintf('export_tikz_panels: matlab2tikz not found, TikZ export of %s skipped\n', base);
    return
end
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

tl = findobj(fig, 'Type', 'tiledlayout');
if isempty(tl)
    axs = findobj(fig, 'Type', 'axes');
else
    axs = findobj(tl(1).Children, 'flat', 'Type', 'axes');
    tiles = arrayfun(@(a) a.Layout.Tile, axs);
    [~, order] = sort(tiles);
    axs = axs(order);
end

for k = 1:numel(axs)
    ax = axs(k);
    nf = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 560 420]);
    if ~isempty(ax.Legend) && isvalid(ax.Legend)
        c = copyobj([ax.Legend, ax], nf);
        nax = c(2);
    else
        nax = copyobj(ax, nf);
    end
    set(nax, 'Units', 'normalized', 'Position', [0.13 0.13 0.8 0.78]);
    nax.XLim = ax.XLim;  nax.YLim = ax.YLim;          % freeze limits

    % xline / yline -> ordinary lines (not supported by matlab2tikz)
    cls = findobj(nax, 'Type', 'ConstantLine');
    for i = 1:numel(cls)
        cl = cls(i);
        if strcmpi(cl.InterceptAxis, 'x')
            xx = [cl.Value cl.Value];  yy = nax.YLim;
        else
            xx = nax.XLim;  yy = [cl.Value cl.Value];
        end
        line(nax, xx, yy, 'Color', cl.Color, 'LineStyle', cl.LineStyle, ...
            'LineWidth', cl.LineWidth, 'HandleVisibility', 'off');
        delete(cl);
    end

    % legend-only dummy lines (all NaN) would be deleted by cleanfigure
    dummies = findobj(nax, 'Type', 'line');
    dummies = dummies(arrayfun(@(l) all(isnan(l.YData)), dummies));
    set(dummies, 'XData', [], 'YData', []);
    lgTitle = '';
    if ~isempty(nax.Legend) && isvalid(nax.Legend)
        lgTitle = nax.Legend.Title.String;
    end

    figure(nf);
    cleanfigure('targetResolution', 300);
    % give the dummies one point outside the (frozen) axes: clipped, but kept
    xo = nax.XLim(1);
    if strcmp(nax.YScale, 'log'), yo = nax.YLim(1) / 100; else, yo = nax.YLim(1) - diff(nax.YLim); end
    set(dummies, 'XData', xo, 'YData', yo);

    file = fullfile(outFolder, sprintf('%s_%d.tex', base, k));
    matlab2tikz(file, 'figurehandle', nf, 'width', '\figurewidth', 'height', '\figureheight', ...
        'showInfo', false, 'showWarnings', false, 'checkForUpdates', false);
    close(nf);

    % legend title -> first legend row (matlab2tikz drops legend titles)
    if ~isempty(lgTitle)
        s = fileread(file);
        i = strfind(s, '\addplot');
        if ~isempty(i)
            s = [s(1:i(1)-1) sprintf('\\addlegendimage{empty legend}\n\\addlegendentry{%s}\n', ...
                lgTitle) s(i(1):end)];
            fid = fopen(file, 'w');  fwrite(fid, s);  fclose(fid);
        end
    end
end
end
