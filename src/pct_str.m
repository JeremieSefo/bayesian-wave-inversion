function s = pct_str(v)
%PCT_STR Format a percentage: one decimal below 10%, none above.
if v < 10, s = sprintf('%.1f', v); else, s = sprintf('%.0f', v); end
end
