function h = plotInstance(params, inst, opts)
%PLOTINSTANCE  Visualize a rendezvous instance (top and 3-D views, side by side).
%
%   h = PLOTINSTANCE(params, inst) draws the region of interest, the denied
%   region, the agents coloured by type, and the weighted centroid in a
%   single figure: a 2-D top view and a 3-D view side by side, sharing one
%   legend. H is the figure handle.
%
%   For a 2-D scenario (params.sim.dim == 2) only the top view is shown.
%   Agent-type labels are taken from PARAMS.types.
%
%   h = PLOTINSTANCE(..., Name, Value) options:
%     'Save'       (false)  write the figure to a file
%     'File'       ('RendezvousDeniedInstance')  output name (no extension)
%     'Format'     ('png')  'png' | 'jpg' | 'tiff' | 'pdf' | 'eps'
%     'Resolution' (300)    dots per inch for raster formats
%   Vector formats ('pdf','eps') are written with vector content; raster
%   formats use the given resolution. The whole side-by-side layout is saved.
%
%   See also GENERATEINSTANCE, DEFAULTPARAMS.

    arguments
        params (1,1) struct
        inst   (1,1) struct
        opts.Save       (1,1) logical = false
        opts.File       (1,:) char    = 'RendezvousDeniedInstance'
        opts.Format     (1,:) char    = 'png'
        opts.Resolution (1,1) double  = 300
    end

    dim = params.sim.dim;
    W   = params.arena.width;
    H   = params.arena.height;
    K   = numel(params.types);

    % ---- colours ------------------------------------------------------
    typeColor     = lines(max(K, 2));             % one distinct colour per type
    roiColor      = [0.9290 0.6940 0.1250];       % orange
    deniedColor   = [0.4940 0.1840 0.5560];       % purple
    centroidColor = [0.8500 0.3250 0.0980];       % red

    roi = inst.roi.vertices;  kr = inst.roi.hull;
    hasDenied = params.sim.denied && ~isempty(inst.denied.vertices);
    if hasDenied
        dR = inst.denied.vertices;  kd = inst.denied.hull;
    end

    % ---- figure + layout ---------------------------------------------
    h = figure('Name', 'Rendezvous Instance', 'Color', 'w');
    nTiles = 1 + (dim == 3);
    t = tiledlayout(h, 1, nTiles, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, 'Rendezvous Instance');

    % ============================================= top view (2-D)
    ax = nexttile(t); hold(ax, 'on'); box(ax, 'on');
    drawRegions(ax, false);
    drawAgents(ax, false);
    axis(ax, 'equal'); xlim(ax, [0 W]); ylim(ax, [0 W]);
    xlabel(ax, 'x (m)'); ylabel(ax, 'y (m)');
    title(ax, 'Top View');
    hold(ax, 'off');

    % ============================================= 3-D view
    if dim == 3
        ax3 = nexttile(t); hold(ax3, 'on'); box(ax3, 'on'); grid(ax3, 'on');
        drawRegions(ax3, true);
        drawAgents(ax3, true);
        xlim(ax3, [0 W]); ylim(ax3, [0 W]); zlim(ax3, [0 H]);
        xlabel(ax3, 'x (m)'); ylabel(ax3, 'y (m)'); zlabel(ax3, 'z (m)');
        view(ax3, [-45 25]);
        title(ax3, '3-D View');
        hold(ax3, 'off');
        legendAx = ax3;
    else
        legendAx = ax;
    end

    % ---- one shared legend for the whole layout ----------------------
    lgd = legend(legendAx, 'Location', 'bestoutside');
    lgd.Layout.Tile = 'east';

    % ---- save the figure -------------------------------------------------
    if opts.Save
        fname = [opts.File '.' opts.Format];
        if ismember(lower(opts.Format), {'pdf', 'eps'})
            exportgraphics(h, fname, 'ContentType', 'vector');
        else
            exportgraphics(h, fname, 'Resolution', opts.Resolution);
        end
        fprintf('Instance figure saved as %s\n', fname);
    end

    % ================================================ nested helpers
    function drawRegions(axh, is3D)
        if is3D
            patch(axh, 'XData', roi(kr,1), 'YData', roi(kr,2), 'ZData', zeros(numel(kr),1), ...
                'FaceColor', roiColor, 'FaceAlpha', 0.08, 'EdgeColor', roiColor, ...
                'LineWidth', 1.5, 'DisplayName', 'Region of interest');
            if hasDenied
                patch(axh, 'XData', dR(kd,1), 'YData', dR(kd,2), 'ZData', zeros(numel(kd),1), ...
                    'FaceColor', deniedColor, 'FaceAlpha', 0.12, 'EdgeColor', deniedColor, ...
                    'LineWidth', 1.5, 'DisplayName', 'Denied region');
            end
        else
            patch(axh, 'XData', roi(kr,1), 'YData', roi(kr,2), ...
                'FaceColor', roiColor, 'FaceAlpha', 0.08, 'EdgeColor', roiColor, ...
                'LineWidth', 1.5, 'DisplayName', 'Region of interest');
            if hasDenied
                patch(axh, 'XData', dR(kd,1), 'YData', dR(kd,2), ...
                    'FaceColor', deniedColor, 'FaceAlpha', 0.12, 'EdgeColor', deniedColor, ...
                    'LineWidth', 1.5, 'DisplayName', 'Denied region');
            end
        end
    end

    function drawAgents(axh, is3D)
        for k = 1:K
            idx = inst.agentType == k;
            if is3D
                scatter3(axh, inst.pos(idx,1), inst.pos(idx,2), inst.pos(idx,3), 55, ...
                    'Marker', '>', 'MarkerFaceColor', typeColor(k,:), ...
                    'MarkerEdgeColor', 'k', 'DisplayName', params.types(k).name);
            else
                scatter(axh, inst.pos(idx,1), inst.pos(idx,2), 55, ...
                    'Marker', '>', 'MarkerFaceColor', typeColor(k,:), ...
                    'MarkerEdgeColor', 'k', 'DisplayName', params.types(k).name);
            end
        end
        if is3D
            plot3(axh, inst.centroid(1), inst.centroid(2), inst.centroid(3), 'o', ...
                'MarkerSize', 10, 'MarkerFaceColor', centroidColor, ...
                'MarkerEdgeColor', 'k', 'DisplayName', 'Weighted centroid');
        else
            plot(axh, inst.centroid(1), inst.centroid(2), 'o', ...
                'MarkerSize', 10, 'MarkerFaceColor', centroidColor, ...
                'MarkerEdgeColor', 'k', 'DisplayName', 'Weighted centroid');
        end
    end
end
