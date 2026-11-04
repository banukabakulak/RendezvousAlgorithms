function h = plotConnectivity(params, results, method, opts)
%PLOTCONNECTIVITY  Windowed joint connectivity (lambda_2) over simulation time.
%
%   h = PLOTCONNECTIVITY(params, results) plots the algebraic connectivity
%   lambda_2 (Fiedler value) of the communication graph aggregated over each
%   sliding window of length T_w = jointWindow/freqInner seconds, against time,
%   using results.conn(r).windowLambda2 from RUNCMAS / RUNDMAS. Each window is
%   placed at the time of its right (ending) edge, so a window at time t
%   summarizes the preceding T_w seconds; lambda_2 > 0 means that window's union
%   graph is connected, so a curve staying above zero is empirical joint
%   connectivity. The convergence threshold (lambda_2 = 0), round boundaries,
%   and the per-run minimum are marked.
%
%   h = PLOTCONNECTIVITY(params, results, method) labels the figure with METHOD
%   ('C-MAS' or 'D-MAS'); when omitted it is inferred from RESULTS.
%
%   Name-Value options (as in PLOTINSTANCE):
%     'Save' (false), 'File' ('JointConnectivity'), 'Format' ('png'),
%     'Resolution' (300).
%
%   See also JOINTCONNECTIVITY, RUNDMAS, RUNCMAS, PLOTPERFORMANCE.

    arguments
        params  (1,1) struct
        results (1,1) struct
        method  (1,:) char = ''
        opts.Save       (1,1) logical = false
        opts.File       (1,:) char    = 'JointConnectivity'
        opts.Format     (1,:) char    = 'png'
        opts.Resolution (1,1) double  = 300
    end

    if isempty(method)
        if isfield(results, 'comm'), method = 'D-MAS'; else, method = 'C-MAS'; end
    end
    if ~isfield(results, 'conn') || isempty(results.conn) || ...
       all(arrayfun(@(c) isempty(c.windowLambda2), results.conn))
        error('plotConnectivity:noData', ...
              ['No windowed connectivity data in results.conn. Run with ' ...
               'params.controller.logGraph = true.']);
    end

    nR    = numel(results.traj);
    dt    = 1 / params.controller.freqInner;
    nStep = size(results.traj{1}.P, 1);
    roundDur = nStep * dt;
    Tw    = params.controller.jointWindow * dt;     % window length, s

    % ---- assemble lambda_2(t) across rounds (window placed at its end time) ----
    time = [];  lam = [];  bnds = [];
    for r = 1:nR
        lw  = results.conn(r).windowLambda2(:);
        win = results.conn(r).window;
        if isempty(lw), continue; end
        ends = (win:win+numel(lw)-1).';             % ending step index per window
        tt   = (r-1) * roundDur + ends * dt;        % global time of the window end
        time = [time; tt];      %#ok<AGROW>
        lam  = [lam;  lw];      %#ok<AGROW>
        bnds = [bnds; r * roundDur]; %#ok<AGROW>
    end

    minLam = min(lam);

    % ---- figure ----
    h  = figure('Name', sprintf('%s Joint Connectivity', method), 'Color', 'w', ...
                'Position', [60 60 900 420]);
    ax = axes(h); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');

    lineCol = [0 0.447 0.741];
    plot(ax, time, lam, '-', 'Color', lineCol, 'LineWidth', 1.5);
    yline(ax, 0, '--k', '\lambda_2 = 0 (disconnected)', 'LineWidth', 1.2, ...
        'FontSize', 10, 'LabelHorizontalAlignment', 'left');

    % round boundaries (skip the final end-of-run boundary)
    for r = 1:nR-1
        xline(ax, r * roundDur, '--', sprintf('Round %d | %d', r, r+1), ...
            'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'FontSize', 9, ...
            'LabelOrientation', 'horizontal');
    end

    xlabel(ax, 'Simulation time (s)');
    ylabel(ax, 'Window algebraic connectivity \lambda_2');
    title(ax, sprintf('%s Joint Connectivity (window T_w = %.0f s)', method, Tw));
    xlim(ax, [0 max(time)]);
    ymax = max(lam) * 1.1; if ymax <= 0, ymax = 1; end
    ylim(ax, [min(0, minLam) ymax]);

    text(ax, 0.02, 0.96, sprintf('min \\lambda_2 = %.3f  (%s)', minLam, ...
        ternary(minLam > 0, 'jointly connected', 'NOT connected')), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, ...
        'BackgroundColor', [0.9 0.95 1], 'EdgeColor', 'k', 'Margin', 3);

    % ---- optional save ----
    if opts.Save
        fname = [opts.File '_' method '.' opts.Format];
        if ismember(lower(opts.Format), {'pdf', 'eps'})
            exportgraphics(h, fname, 'ContentType', 'vector');
        else
            exportgraphics(h, fname, 'Resolution', opts.Resolution);
        end
        fprintf('Connectivity figure saved as %s\n', fname);
    end
end

% =====================================================================
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
