function h = plotPerformance(params, results, method)
%PLOTPERFORMANCE  Controller performance: position convergence and energy.
%
%   h = PLOTPERFORMANCE(params, results) draws a 1x2 figure from the RESULTS
%   of RUNCMAS or RUNDMAS:
%     left   distance of each agent to the round rendezvous point r*_s over
%            time, with the convergence tolerance and round boundaries marked;
%     right  cumulative energy (mAh) of each agent over time.
%
%   h = PLOTPERFORMANCE(params, results, method) labels the figure and title
%   with METHOD ('C-MAS' or 'D-MAS'). When METHOD is omitted it is inferred
%   from RESULTS (the distributed run carries a per-round .comm field).
%
%   Agents are tracked by their persistent id, so additions and removals
%   appear as series that begin or end mid-run (dashed when not present from
%   the start). Colours encode agent type (Type-1 blue, Type-2 red).
%
%   See also RUNCMAS, RUNDMAS, COMPUTEMETRICS, COMPUTEMETRICSDMAS.

    arguments
        params  (1,1) struct
        results (1,1) struct
        method  (1,:) char = ''
    end

    % method label: explicit argument, else inferred from the results struct
    if isempty(method)
        if isfield(results, 'comm'), method = 'D-MAS'; else, method = 'C-MAS'; end
    end

    nR    = numel(results.traj);
    nStep = size(results.traj{1}.P, 1);
    dt    = 1 / params.controller.freqInner;
    roundDur = nStep * dt;
    tol   = params.controller.epsilon;          % convergence tolerance, m

    % ---- global timeline ----
    time = zeros(1, nR * nStep);
    for r = 1:nR
        time((r-1)*nStep + (1:nStep)) = (0:nStep-1) * dt + (r-1) * roundDur;
    end
    Ttot = nR * nStep;

    % ---- collect persistent agents and their types ----
    allID = [];
    for r = 1:nR, allID = union(allID, results.traj{r}.id); end
    allID = allID(:);
    nA    = numel(allID);
    gtype = zeros(nA, 1);
    present1 = false(nA, 1);                     % present in round 1?
    for a = 1:nA
        for r = 1:nR
            loc = find(results.traj{r}.id == allID(a), 1);
            if ~isempty(loc)
                gtype(a) = results.traj{r}.type(loc);
                if r == 1, present1(a) = true; end
                break;
            end
        end
    end

    % ---- assemble distance and cumulative-energy series per agent ----
    dist = nan(nA, Ttot);
    ener = nan(nA, Ttot);
    for r = 1:nR
        P   = results.traj{r}.P;                 % nStep x Nr x dim
        E   = results.traj{r}.energy;            % nStep x Nr
        ids = results.traj{r}.id;
        rStar = results.rendezvous(r, :);
        cols = (r-1)*nStep + (1:nStep);
        for i = 1:numel(ids)
            a = find(allID == ids(i), 1);
            Pi = reshape(P(:, i, :), nStep, []);         % nStep x dim
            dist(a, cols) = vecnorm(Pi - rStar, 2, 2).';
            ener(a, cols) = E(:, i).';
        end
    end

    % ---- colours by type (Type-1 blue, Type-2 red) ----
    idx1 = find(gtype == 1);  idx2 = find(gtype == 2);
    blueT1 = [0 0.447 0.741];      % Type-1 (scout)
    redT2  = [0.85 0.325 0.098];   % Type-2 (carrier/relay)
    color = zeros(nA, 3);
    color(idx1, :) = repmat(blueT1, numel(idx1), 1);
    color(idx2, :) = repmat(redT2,  numel(idx2), 1);

    % ---- figure ----
    h = figure('Name', sprintf('%s Performance', method), 'Color', 'w', 'Position', [60 60 1500 560]);
    t = tiledlayout(h, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('%s Controller Performance', method), 'FontWeight', 'bold');

    % ===================== left: position convergence =====================
    ax1 = nexttile(t); hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');
    yTop = max(dist(:), [], 'omitnan') * 1.1;
    % sustained-convergence interval of each round shaded green (vertical band
    % from the convergence time T_c to the round end), behind the series
    for r = 1:nR
        tc = results.consensus(r).time;
        if ~isnan(tc)
            bStart = (r-1) * roundDur + tc;
            bEnd   = r * roundDur;
            patch(ax1, [bStart bEnd bEnd bStart], [0 0 yTop yTop], [0.20 0.70 0.30], ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    plotSeries(ax1, time, dist, color, present1);
    yline(ax1, tol, '--k', sprintf('\\epsilon = %g m', tol), ...
        'LineWidth', 1.3, 'FontSize', 10, 'LabelHorizontalAlignment', 'left');
    drawRoundMarkers(ax1, results, nR, roundDur, time, dist);
    xlabel(ax1, 'Simulation time (s)'); ylabel(ax1, 'Distance to r^*_s (m)');
    title(ax1, 'Position Convergence');
    xlim(ax1, [0 max(time)]); ylim(ax1, [0 yTop]);
    addTypeLegend(ax1, idx1, idx2, present1, color, 'northeast');

    % ===================== right: cumulative energy =======================
    ax2 = nexttile(t); hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
    plotSeries(ax2, time, ener, color, present1);
    drawRoundMarkers(ax2, results, nR, roundDur, time, ener);
    xlabel(ax2, 'Simulation time (s)'); ylabel(ax2, 'Cumulative energy (mAh)');
    title(ax2, 'Energy Consumption');
    xlim(ax2, [0 max(time)]); ylim(ax2, [0 max(ener(:), [], 'omitnan') * 1.1]);
    addTypeLegend(ax2, idx1, idx2, present1, color, 'northwest');

    % average final energy per type (annotation)
    finalE = ener(:, end);
    if ~isempty(idx1)
        text(ax2, max(time)*0.75, max(ener(:),[],'omitnan')*0.95, ...
            sprintf('Type-1 avg: %.1f mAh', mean(finalE(idx1), 'omitnan')), ...
            'FontSize', 9, 'BackgroundColor', [0.85 0.92 1], 'EdgeColor', 'k', 'Margin', 2);
    end
    if ~isempty(idx2)
        text(ax2, max(time)*0.75, max(ener(:),[],'omitnan')*0.86, ...
            sprintf('Type-2 avg: %.1f mAh', mean(finalE(idx2), 'omitnan')), ...
            'FontSize', 9, 'BackgroundColor', [1 0.85 0.85], 'EdgeColor', 'k', 'Margin', 2);
    end
end

% =====================================================================
function plotSeries(ax, time, Y, color, present1)
    for a = 1:size(Y, 1)
        ls = '-'; if ~present1(a), ls = ':'; end
        plot(ax, time, Y(a, :), 'Color', color(a, :), 'LineWidth', 1.4, 'LineStyle', ls);
    end
end

% =====================================================================
function drawRoundMarkers(ax, results, nR, roundDur, time, Y)
% Vertical line and label at each inter-round boundary (labelled by event),
% plus a "Round r" tag centred in each round.
    ymax = max(Y(:), [], 'omitnan');
    for r = 1:nR-1
        tb = r * roundDur;
        cur = results.traj{r}.id;  nxt = results.traj{r+1}.id;
        added   = ~isempty(setdiff(nxt, cur));
        removed = ~isempty(setdiff(cur, nxt));
        if added && removed,   lbl = 'Add + Remove';
        elseif added,          lbl = 'Agent Addition';
        elseif removed,        lbl = 'Agent Removal';
        else,                  lbl = '';
        end
        xline(ax, tb, '--', lbl, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
            'FontSize', 10, 'LabelOrientation', 'horizontal');
    end
    for r = 1:nR
        text(ax, (r-0.5)*roundDur, ymax*0.92, sprintf('Round %d', r), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
            'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 10);
    end
end

% =====================================================================
function addTypeLegend(ax, idx1, idx2, present1, color, loc)
    hs = []; lb = {};
    if ~isempty(idx1)
        hs(end+1) = plot(ax, NaN, NaN, '-', 'Color', color(idx1(1), :), 'LineWidth', 1.6);
        lb{end+1} = 'Type-1 (scout)';
    end
    if ~isempty(idx2)
        hs(end+1) = plot(ax, NaN, NaN, '-', 'Color', color(idx2(1), :), 'LineWidth', 1.6);
        lb{end+1} = 'Type-2 (carrier/relay)';
    end
    if any(~present1)
        hs(end+1) = plot(ax, NaN, NaN, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.6);
        lb{end+1} = 'Added mid-run';
    end
    if ~isempty(hs), legend(ax, hs, lb, 'Location', loc, 'FontSize', 9); end
end
