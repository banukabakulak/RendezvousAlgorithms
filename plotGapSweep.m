function h = plotGapSweep(out, whichGap, opts)
%PLOTGAPSWEEP  D-MAS/C-MAS gap versus topology sparsity (arena size and diameter).
%
%   h = PLOTGAPSWEEP(out) plots the D-MAS optimality gap from EXPERIMENTGAPSWEEP
%   against the arena side L (left) and the initial graph diameter (right).
%
%   h = PLOTGAPSWEEP(out, whichGap) selects the gap metric: 'obj' (default, the
%   CRM-objective optimality gap) or 'energy' (the mission-energy gap, available
%   only if EXPERIMENTGAPSWEEP was run with doEnergy = true).
%
%   Name-Value options (as in PLOTINSTANCE): 'Save' (false), 'File' ('GapSweep'),
%   'Format' ('png'), 'Resolution' (300).
%
%   See also EXPERIMENTGAPSWEEP.

arguments
    out      (1,1) struct
    whichGap (1,:) char {mustBeMember(whichGap,{'obj','ptDev','energy'})} = 'obj'
    opts.Layout     (1,:) char {mustBeMember(opts.Layout,{'side','stacked'})} = 'side'
    opts.Save       (1,1) logical = false
    opts.File       (1,:) char    = 'GapSweep'
    opts.Format     (1,:) char    = 'png'
    opts.Resolution (1,1) double  = 300
end

switch whichGap
    case 'ptDev'
        gap = out.ptDev;      gname = 'Point deviation (% of arena side)';
    case 'obj'
        gap = out.objGap;     gname = 'Optimality gap (CRM objective, %)';
    case 'energy'
        gap = out.energyGap;  gname = 'Mission-energy gap (%)';
        if all(isnan(gap))
            error('plotGapSweep:noEnergy', ...
                'No energy gap in out; run experimentGapSweep(..., true).');
        end
end

col = [0 0.447 0.741];
if strcmp(opts.Layout, 'stacked')       % single-column: two panels stacked
    pos = [60 60 460 620];  rc = [2 1];
else                                      % two-column: side by side
    pos = [60 60 1000 400]; rc = [1 2];
end
h = figure('Name', 'D-MAS Optimality Gap vs Sparsity', 'Color', 'w', 'Position', pos);
t = tiledlayout(h, rc(1), rc(2), 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(t); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, out.L, gap, '-o', 'Color', col, 'MarkerFaceColor', col, 'LineWidth', 1.6);
xlabel(ax1, 'Arena side L (m)'); ylabel(ax1, gname);
title(ax1, 'Gap vs arena size');

ax2 = nexttile(t); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, out.diam, gap, '-o', 'Color', col, 'MarkerFaceColor', col, 'LineWidth', 1.6);
xlabel(ax2, 'Initial graph diameter (hops)'); ylabel(ax2, gname);
title(ax2, 'Gap vs diameter');

title(t, 'D-MAS Optimality Gap versus Topology Sparsity', 'FontWeight', 'bold');

if opts.Save
    fname = [opts.File '.' opts.Format];
    if ismember(lower(opts.Format), {'pdf','eps'})
        exportgraphics(h, fname, 'ContentType', 'vector');
    else
        exportgraphics(h, fname, 'Resolution', opts.Resolution);
    end
    fprintf('Gap-sweep figure saved as %s\n', fname);
end
end
