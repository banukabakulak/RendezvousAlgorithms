function h = plotConnectivityCompare(params, resultsC, resultsD, opts)
%PLOTCONNECTIVITYCOMPARE  Compare C-MAS and D-MAS joint connectivity on one axis.
%
%   h = PLOTCONNECTIVITYCOMPARE(params, resultsC, resultsD) overlays the
%   windowed algebraic connectivity lambda_2(t) of the C-MAS run (RESULTSC, from
%   RUNCMAS) and the D-MAS run (RESULTSD, from RUNDMAS) on a single time axis,
%   with the lambda_2 = 0 threshold, round boundaries, and each run's minimum
%   marked. lambda_2 is computed over sliding windows of length
%   T_w = jointWindow/freqInner seconds (see JOINTCONNECTIVITY); a curve above
%   zero everywhere is empirical joint connectivity.
%
%   Name-Value options (as in PLOTINSTANCE):
%     'Save' (false), 'File' ('JointConnectivity_compare'), 'Format' ('png'),
%     'Resolution' (300).
%
%   See also PLOTCONNECTIVITY, JOINTCONNECTIVITY, RUNCMAS, RUNDMAS.

arguments
    params   (1,1) struct
    resultsC (1,1) struct
    resultsD (1,1) struct
    opts.Save       (1,1) logical = false
    opts.File       (1,:) char    = 'JointConnectivity_compare'
    opts.Format     (1,:) char    = 'png'
    opts.Resolution (1,1) double  = 300
end

dt = 1 / params.controller.freqInner;
Tw = params.controller.jointWindow * dt;            % window length, s

[tC, lamC, roundDur, nR] = buildSeries(params, resultsC, dt);
[tD, lamD, ~,        ~ ] = buildSeries(params, resultsD, dt);

colC = [0 0.447 0.741];        % C-MAS blue
colD = [0.85 0.325 0.098];     % D-MAS red

h  = figure('Name', 'Joint Connectivity: C-MAS vs D-MAS', 'Color', 'w', ...
    'Position', [60 60 900 420]);
ax = axes(h); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');

pC = plot(ax, tC, lamC, '-',  'Color', colC, 'LineWidth', 1.6);
pD = plot(ax, tD, lamD, '-',  'Color', colD, 'LineWidth', 1.6);
yline(ax, 0, '--k', '\lambda_2 = 0 (disconnected)', 'LineWidth', 1.2, ...
    'FontSize', 10, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');

for r = 1:nR-1
    xline(ax, r * roundDur, '--', sprintf('Round %d | %d', r, r+1), ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'FontSize', 9, ...
        'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
end

xlabel(ax, 'Simulation time (s)');
ylabel(ax, 'Window algebraic connectivity \lambda_2');
title(ax, sprintf('Joint Connectivity (window T_w = %.0f s)', Tw));

tMax = max([tC; tD]);  xlim(ax, [0 tMax]);
ymax = max([lamC; lamD]) * 1.1; if ymax <= 0, ymax = 1; end
ymin = min([0; lamC; lamD]);
ylim(ax, [ymin ymax]);

legend(ax, [pC pD], ...
    {sprintf('C-MAS (min \\lambda_2 = %.3f)', min(lamC)), ...
    sprintf('D-MAS (min \\lambda_2 = %.3f)', min(lamD))}, ...
    'Location', 'southeast', 'FontSize', 10);

if opts.Save
    fname = [opts.File '.' opts.Format];
    if ismember(lower(opts.Format), {'pdf', 'eps'})
        exportgraphics(h, fname, 'ContentType', 'vector');
    else
        exportgraphics(h, fname, 'Resolution', opts.Resolution);
    end
    fprintf('Connectivity comparison figure saved as %s\n', fname);
end
end

% =====================================================================
function [time, lam, roundDur, nR] = buildSeries(params, results, dt)
%BUILDSERIES  Assemble lambda_2(t) across rounds (window placed at its end time).
if ~isfield(results, 'conn') || isempty(results.conn)
    error('plotConnectivityCompare:noData', ...
        ['No windowed connectivity data in results.conn. Run with ' ...
        'params.controller.logGraph = true.']);
end
nR       = numel(results.traj);
nStep    = size(results.traj{1}.P, 1);
roundDur = nStep * dt;
time = [];  lam = [];
for r = 1:nR
    lw  = results.conn(r).windowLambda2(:);
    win = results.conn(r).window;
    if isempty(lw), continue; end
    ends = (win:win+numel(lw)-1).';
    time = [time; (r-1) * roundDur + ends * dt];  %#ok<AGROW>
    lam  = [lam;  lw];                             %#ok<AGROW>
end
end
