function out = experimentComparisonFeasibility(nInst, baseParams)
%EXPERIMENTCOMPARISONFEASIBILITY  Feasibility and energy-optimality of the
%   rendezvous POINT chosen by each method, with the denied region ACTIVE.
%
%   out = EXPERIMENTCOMPARISONFEASIBILITY(nInst, baseParams) complements
%   EXPERIMENTCOMPARISON by exposing what the baselines do NOT do: they only
%   select a rendezvous point, without ensuring it is feasible and without
%   accounting for energy. With the denied region enabled (pinned to the weighted
%   centroid, as in the main scenario), for each method we report, over nInst
%   instances:
%
%     feasible (%)  fraction of instances whose chosen point lies in
%                   R \ D (inside the RoI and outside the denied region)
%     Delta f (%)   CRM-objective gap of the point to the energy-weighted
%                   optimum r* (the C-MAS / H-CRM-D point):
%                   100 * (f(r) - f(r*)) / f(r*),  f(r)=sum_i c_i||p_i-r||^2
%
%   The proposed methods solve the constrained, energy-weighted problem, so they
%   are feasible by construction (C-MAS exactly optimal, D-MAS near-optimal). The
%   baselines only place a point:
%     CA  unweighted centroid       (ignores energy; no feasibility handling)
%     RC  RoI geometric center      (ignores agents and energy; no feasibility)
%     CC  minimax circumcenter      (optimizes fairness, not energy; no feasibility)
%
%   Defaults: nInst = 300, baseParams = defaultParams().
%
%   OUT has fields method, feasible (%), gap (mean Delta f, %), gapFeas (mean
%   Delta f over feasible instances only), and OUT.raw.
%
%   See also EXPERIMENTCOMPARISON, SOLVERENDEZVOUSCENTRAL, SOLVERENDEZVOUSDISTRIBUTED.

    arguments
        nInst      (1,1) double = 300
        baseParams (1,1) struct = defaultParams()
    end

    p = baseParams;
    p.sim.N                   = 30;
    p.sim.denied              = true;      % expose feasibility: denied region active
    p.sim.forceDeniedCentroid = true;      % pinned to the weighted centroid (main scenario)
    p.controller.logGraph     = false;

    method = {'C-MAS','D-MAS','CA','RC','CC'};
    nM = numel(method);
    feas = nan(nInst,nM);  gap = nan(nInst,nM);

    for s = 1:nInst
        p.sim.seed = 15000 + s;
        inst = generateInstance(p);
        P = inst.pos;  w = inst.weight;  dim = size(P,2);
        f = @(r) sum(w .* sum((P - r).^2, 2));

        rOpt = solveRendezvousCentral(p, inst);               % C-MAS feasible optimum
        cg   = buildCommGraph(p, inst);
        rD   = solveRendezvousDistributed(p, inst, cg);        % D-MAS point

        pts = {rOpt(:).', rD(:).', ...
               mean(P,1), ...                                  % CA
               roiCenter(inst, P, dim), ...                    % RC
               circumCenter(P, dim)};                          % CC

        fOpt = max(f(rOpt), eps);
        for mi = 1:nM
            r = pts{mi};
            feas(s,mi) = isFeasible(r, inst);
            gap(s,mi)  = 100 * (f(r) - fOpt) / fOpt;
        end
    end

    out.method  = method(:);
    out.feasible = 100 * mean(feas, 1, 'omitnan').';
    out.gap      = mean(gap, 1, 'omitnan').';
    gf = nan(nM,1);
    for mi = 1:nM
        m = feas(:,mi) == 1;
        if any(m), gf(mi) = mean(gap(m,mi), 'omitnan'); end
    end
    out.gapFeas = gf;
    out.raw = struct('feas',feas,'gap',gap);

    printTable(out);
    printLatex(out);
end

% =====================================================================
function tf = isFeasible(r, inst)
%ISFEASIBLE  Point in R \ D: inside the RoI polygon and outside the denied one.
    x = r(1);  y = r(2);
    V = inst.roi.vertices;
    inRoI = inpolygon(x, y, V(:,1), V(:,2));
    inDen = false;
    if isfield(inst,'denied') && ~isempty(inst.denied.vertices)
        D = inst.denied.vertices;
        inDen = inpolygon(x, y, D(:,1), D(:,2));
    end
    tf = inRoI && ~inDen;
end

% ---------------------------------------------------------------------
function r = roiCenter(inst, P, dim)
    cxy = mean(inst.roi.vertices, 1);
    if dim >= 3, r = [cxy, mean(P(:,3))]; else, r = cxy; end
end

% ---------------------------------------------------------------------
function r = circumCenter(P, dim)
    c0 = mean(P(:,1:2), 1);
    fc = @(c) max(vecnorm(P(:,1:2) - c, 2, 2));
    c  = fminsearch(fc, c0, optimset('Display','off','TolX',1e-4,'TolFun',1e-4));
    if dim >= 3, r = [c, mean(P(:,3))]; else, r = c; end
end

% ---------------------------------------------------------------------
function printTable(out)
    fprintf('\n=== Feasibility & energy-optimality of the rendezvous point ===\n');
    fprintf('(denied region active, %d instances)\n', numel(out.raw.feas(:,1)));
    fprintf('%-8s %-12s %-14s %-14s\n', 'Method','feasible(%)','Df all(%)','Df feas-only(%)');
    for i = 1:numel(out.method)
        fprintf('%-8s %-12.1f %-14.2f %-14.2f\n', ...
            out.method{i}, out.feasible(i), out.gap(i), out.gapFeas(i));
    end
    fprintf('===============================================================\n');
end

% ---------------------------------------------------------------------
function printLatex(out)
    fprintf('\n%% --- LaTeX rows for tab:ComparisonFeasibility ---\n');
    for i = 1:numel(out.method)
        fprintf('\\texttt{%s} & $%.1f$ & $%.2f$ \\\\\n', ...
            out.method{i}, out.feasible(i), out.gap(i));
    end
    fprintf('%% ------------------------------------------------\n');
end
