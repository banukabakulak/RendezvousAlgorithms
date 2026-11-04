function animateRendezvous(params, results, opts)
%ANIMATERENDEZVOUS  Animate MAS rendezvous (3-D, 2-D, zoomed-2-D) per round.
%
%   ANIMATERENDEZVOUS(params, results) animates the trajectories from RUNCMAS
%   in a (nRounds x 3) grid: 3-D view, top view, and a zoomed top view around
%   the rendezvous point, one row per round.
%
%   ANIMATERENDEZVOUS(..., Name, Value) options:
%     'Record'     (false)  write an .avi/.mp4 video of the animation
%     'File'       ('MAS_rendezvous')  output file name (no extension)
%     'FrameRate'  (30)      video frame rate
%     'ZoomRadius' (1.5)     half-width of the zoomed view, m
%     'PauseTime'  (0.01)    pause between frames when not recording, s
%     'Stride'     (1)       animate every Stride-th control step (speed-up)
%
%   Agents are coloured by type; agents added mid-run are highlighted. The
%   region of interest and denied region are drawn in every panel.
%
%   See also RUNCMAS, PLOTTRAJECTORIES.

    arguments
        params  (1,1) struct
        results (1,1) struct
        opts.Record     (1,1) logical = false
        opts.File       (1,:) char    = 'MAS_rendezvous'
        opts.FrameRate  (1,1) double  = 30
        opts.ZoomRadius (1,1) double  = 1.5
        opts.PauseTime  (1,1) double  = 0.01
        opts.Stride     (1,1) double  = 1
    end

    nR  = numel(results.traj);
    W   = params.arena.width;
    roi = results.instFinal.roi;
    denied = results.instFinal.denied;
    hasDenied = params.sim.denied && ~isempty(denied.vertices);

    roiC = [0.9290 0.6940 0.1250];
    denC = [0.4940 0.1840 0.5560];
    rvC  = [1 0.5 0];
    typeCol = [0 0.447 0.741; 0.85 0.325 0.098];
    newC = [1 0.5 0];
    id1  = results.traj{1}.id;                 % ids present from the start

    v = [];
    if opts.Record
        v = VideoWriter(opts.File);            %#ok<TNMLP>
        v.FrameRate = opts.FrameRate;
        open(v);
    end

    fig = figure('Name', 'MAS Rendezvous Animation', 'Color', 'w', ...
        'Position', [40 40 1200 420*nR]);

    for r = 1:nR
        P   = results.traj{r}.P;               % nStep x Nr x dim
        ids = results.traj{r}.id;
        typ = results.traj{r}.type;
        Nr  = size(P, 2);
        nStep = size(P, 1);
        rv  = results.rendezvous(r, :);

        i3D = (r-1)*3 + 1;  i2D = (r-1)*3 + 2;  iZ = (r-1)*3 + 3;

        % --- set up the three panels for this round ---
        ax3 = subplot(nR, 3, i3D); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
        title(ax3, sprintf('Round %d — 3D', r)); view(ax3,[-45 25]);
        axis(ax3,[0 W 0 W 0 params.arena.height]); axis(ax3,'vis3d');
        xlabel(ax3,'x'); ylabel(ax3,'y'); zlabel(ax3,'z');

        ax2 = subplot(nR, 3, i2D); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
        title(ax2, sprintf('Round %d — Top', r)); axis(ax2,'equal'); axis(ax2,[0 W 0 W]);
        xlabel(ax2,'x'); ylabel(ax2,'y');

        axZ = subplot(nR, 3, iZ); hold(axZ,'on'); grid(axZ,'on'); box(axZ,'on');
        title(axZ, sprintf('Round %d — Zoom', r)); axis(axZ,'equal');
        axis(axZ,[rv(1)-opts.ZoomRadius rv(1)+opts.ZoomRadius ...
                  rv(2)-opts.ZoomRadius rv(2)+opts.ZoomRadius]);
        xlabel(axZ,'x'); ylabel(axZ,'y');

        % region overlays + rendezvous marker in each panel
        drawRegion3(ax3, roi, roiC); drawRegion2(ax2, roi, roiC); drawRegion2(axZ, roi, roiC);
        if hasDenied
            drawRegion3(ax3, denied, denC); drawRegion2(ax2, denied, denC); drawRegion2(axZ, denied, denC);
        end
        plot3(ax3, rv(1), rv(2), rv(3), 'p', 'MarkerSize', 12, 'MarkerFaceColor', rvC, 'MarkerEdgeColor', rvC);
        plot(ax2, rv(1), rv(2), 'p', 'MarkerSize', 12, 'MarkerFaceColor', rvC, 'MarkerEdgeColor', rvC);
        plot(axZ, rv(1), rv(2), 'p', 'MarkerSize', 14, 'MarkerFaceColor', rvC, 'MarkerEdgeColor', rvC);

        % per-agent animated handles (trajectory + current marker)
        tr3 = gobjects(Nr,1); cu3 = gobjects(Nr,1);
        tr2 = gobjects(Nr,1); cu2 = gobjects(Nr,1);
        trZ = gobjects(Nr,1); cuZ = gobjects(Nr,1);
        for i = 1:Nr
            col = typeCol(typ(i), :);
            if ~ismember(ids(i), id1), col = newC; end     % added mid-run
            tr3(i) = plot3(ax3, NaN, NaN, NaN, '.', 'Color', col, 'MarkerSize', 4);
            cu3(i) = plot3(ax3, NaN, NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'MarkerSize', 5);
            tr2(i) = plot(ax2, NaN, NaN, '.', 'Color', col, 'MarkerSize', 4);
            cu2(i) = plot(ax2, NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'MarkerSize', 5);
            trZ(i) = plot(axZ, NaN, NaN, '.', 'Color', col, 'MarkerSize', 4);
            cuZ(i) = plot(axZ, NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'MarkerSize', 6);
        end

        % --- animate ---
        for tk = 1:opts.Stride:nStep
            for i = 1:Nr
                x = P(1:tk, i, 1); y = P(1:tk, i, 2); z = P(1:tk, i, 3);
                set(tr3(i), 'XData', x, 'YData', y, 'ZData', z);
                set(cu3(i), 'XData', x(end), 'YData', y(end), 'ZData', z(end));
                set(tr2(i), 'XData', x, 'YData', y);
                set(cu2(i), 'XData', x(end), 'YData', y(end));
                set(trZ(i), 'XData', x, 'YData', y);
                set(cuZ(i), 'XData', x(end), 'YData', y(end));
            end
            drawnow;
            if opts.Record
                writeVideo(v, getframe(fig));
            elseif opts.PauseTime > 0
                pause(opts.PauseTime);
            end
        end

        % final-position stars
        for i = 1:Nr
            plot3(ax3, P(end,i,1), P(end,i,2), P(end,i,3), '*', 'Color', 'r', 'MarkerSize', 6);
            plot(ax2, P(end,i,1), P(end,i,2), '*', 'Color', 'r', 'MarkerSize', 6);
            plot(axZ, P(end,i,1), P(end,i,2), '*', 'Color', 'r', 'MarkerSize', 7);
        end
    end

    if opts.Record
        close(v);
        fprintf('Video saved as %s\n', opts.File);
    end
end

% =====================================================================
function drawRegion3(ax, reg, color)
    V = reg.vertices; k = reg.hull;
    plot3(ax, V(k,1), V(k,2), zeros(numel(k),1), '-', 'Color', color, 'LineWidth', 1.2);
end
function drawRegion2(ax, reg, color)
    V = reg.vertices; k = reg.hull;
    plot(ax, V(k,1), V(k,2), '-', 'Color', color, 'LineWidth', 1.2);
end
