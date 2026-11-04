function h = plotTrajectories(params, results, method)
%PLOTTRAJECTORIES  Per-round 3-D agent trajectories with rendezvous points.
%
%   h = PLOTTRAJECTORIES(params, results) draws one 3-D panel per round from
%   the RESULTS of RUNCMAS or RUNDMAS: each agent's start (circle), path
%   (type-coloured), and final position (star), with the region of interest,
%   the denied region, and the round rendezvous point. H is the figure handle.
%
%   h = PLOTTRAJECTORIES(params, results, method) labels the figure and title
%   with METHOD ('C-MAS' or 'D-MAS'). When METHOD is omitted it is inferred
%   from RESULTS (the distributed run carries a per-round .comm field).
%
%   See also RUNCMAS, RUNDMAS, PLOTPERFORMANCE, ANIMATERENDEZVOUS.

arguments
    params  (1,1) struct
    results (1,1) struct
    method  (1,:) char = ''
end

% method label: explicit argument, else inferred from the results struct
if isempty(method)
    if isfield(results, 'comm'), method = 'D-MAS'; else, method = 'C-MAS'; end
end

nR  = numel(results.traj);
W   = params.arena.width;
H   = params.arena.height;
roi = results.instFinal.roi;
denied = results.instFinal.denied;

roiC    = [0.9290 0.6940 0.1250];
denC    = [0.4940 0.1840 0.5560];
rvC     = [1 0.5 0];
typeCol = [0 0.447 0.741; 0.85 0.325 0.098];   % type-1, type-2

h = figure('Name', sprintf('%s Trajectories', method), 'Color', 'w', 'Position', [60 60 560*nR 560]);
t = tiledlayout(h, 1, nR, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, sprintf('%s Rendezvous Trajectories', method), 'FontWeight', 'bold');

for r = 1:nR
    ax = nexttile(t); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');

    drawRegion(ax, roi,  roiC);
    if params.sim.denied && ~isempty(denied.vertices)
        drawRegion(ax, denied, denC);
    end

    P    = results.traj{r}.P;              % nStep x Nr x dim
    typ  = results.traj{r}.type;
    Nr   = size(P, 2);
    for i = 1:Nr
        xyz = reshape(P(:, i, :), size(P,1), []);     % nStep x dim
        col = typeCol(typ(i), :);
        plot3(ax, xyz(:,1), xyz(:,2), xyz(:,3), '-', 'Color', col, 'LineWidth', 0.8);
        plot3(ax, xyz(1,1), xyz(1,2), xyz(1,3), 'o', 'Color', 'g', ...
            'MarkerFaceColor', 'g', 'MarkerSize', 5);
        plot3(ax, xyz(end,1), xyz(end,2), xyz(end,3), '*', 'Color', 'r', 'MarkerSize', 6);
    end

    rv = results.rendezvous(r, :);
    plot3(ax, rv(1), rv(2), rv(3), 'p', 'MarkerSize', 14, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', rvC, 'MarkerFaceColor', rvC);

    axis(ax, [0 W 0 W 0 H]);
    view(ax, [-45 25]);
    xlabel(ax, 'x (m)'); ylabel(ax, 'y (m)'); zlabel(ax, 'z (m)');
    title(ax, sprintf('Round %d  (N = %d)', r, Nr));
end
end

% =====================================================================
function drawRegion(ax, reg, color)
V = reg.vertices;  k = reg.hull;
patch(ax, 'XData', V(k,1), 'YData', V(k,2), 'ZData', zeros(numel(k),1), ...
    'FaceColor', color, 'FaceAlpha', 0.08, 'EdgeColor', color, 'LineWidth', 1.2);
end
