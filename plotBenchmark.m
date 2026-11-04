function h = plotBenchmark(out, opts)
%PLOTBENCHMARK  Log-log runtime of CRM-OPT versus general QP solvers.
%
%   h = PLOTBENCHMARK(out) plots the mean wall-clock time (ms) of each solver in
%   OUT (from EXPERIMENTBENCHMARK) against the region size M on a log-log scale.
%   CRM-OPT is essentially flat in M while the interior-point solver grows steeply;
%   the markers confirm all solvers return the same optimum.
%
%   Name-Value options (as in PLOTINSTANCE): 'Save' (false), 'File' ('Benchmark'),
%   'Format' ('png'), 'Resolution' (300).
%
%   See also EXPERIMENTBENCHMARK.

arguments
    out              (1,1) struct
    opts.Save        (1,1) logical = false
    opts.File        (1,:) char    = 'Benchmark'
    opts.Format      (1,:) char    = 'png'
    opts.Resolution  (1,1) double  = 300
end

cCRM = [0 0.447 0.741];      % blue
cAS  = [0.466 0.674 0.188];  % green
cIP  = [0.85 0.325 0.098];   % red
cCom = [0.494 0.184 0.556];  % purple

h = figure('Name','CRM-OPT Runtime Benchmark','Color','w','Position',[80 80 640 460]);
ax = axes(h); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'XScale','log','YScale','log');

loglog(ax, out.M, out.tCRM, '-o', 'Color',cCRM, 'MarkerFaceColor',cCRM, 'LineWidth',1.8);
loglog(ax, out.M, out.tAS,  '-s', 'Color',cAS,  'MarkerFaceColor',cAS,  'LineWidth',1.6);
loglog(ax, out.M, out.tIP,  '-^', 'Color',cIP,  'MarkerFaceColor',cIP,  'LineWidth',1.6);
names = {'\texttt{CRM-OPT}','\texttt{quadprog}-AS','\texttt{quadprog}-IP'};
if isfield(out,'hasCom') && out.hasCom
    loglog(ax, out.M, out.tCom, '-d', 'Color',cCom, 'MarkerFaceColor',cCom, 'LineWidth',1.6);
    names{end+1} = out.comName;
end

xlabel(ax, 'Region size $M$ (facets / constraints)', 'Interpreter','latex');
ylabel(ax, 'Wall-clock time (ms)', 'Interpreter','latex');
title(ax,  '\texttt{CRM-OPT} \ vs.\ general QP solvers', 'Interpreter','latex');
legend(ax, names, 'Interpreter','latex', 'Location','northwest');
xlim(ax, [min(out.M)*0.8, max(out.M)*1.25]);

if opts.Save
    fname = [opts.File '.' opts.Format];
    if ismember(lower(opts.Format), {'pdf','eps'})
        exportgraphics(h, fname, 'ContentType','vector');
    else
        exportgraphics(h, fname, 'Resolution', opts.Resolution);
    end
    fprintf('Benchmark figure saved as %s\n', fname);
end
end
