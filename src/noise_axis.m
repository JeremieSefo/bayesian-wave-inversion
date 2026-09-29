function noise_axis(ax, kEff)
%NOISE_AXIS Label a logarithmic sigma axis with the effective relative noise.
%   Ticks at sigma = 0.01, 0.1, 1 with the effective relative noise
%   100*kEff*sigma (in %) written underneath.
t = [0.01 0.1 1];
base = {'10^{-2}', '10^{-1}', '10^{0}'};
lbl = cell(size(t));
for i = 1:numel(t)
    v = 100 * kEff * t(i);
    pct = pct_str(v);
    lbl{i} = sprintf('%s\\newline(%s%%)', base{i}, pct);
end
set(ax, 'XTick', t, 'XTickLabel', lbl, 'TickLabelInterpreter', 'tex');
end
