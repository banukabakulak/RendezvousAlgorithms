function out = experimentDeniedDensity(coverage, nC, nD, opts)
%EXPERIMENTDENIEDDENSITY  Denied-region density and overlap analysis (Sec. density).
%
%   out = EXPERIMENTDENIEDDENSITY(coverage, nC, nD) sweeps the denied-region
%   density and reports both tables of Sec.~\ref{sec:density}:
%
%   (1) Centralized H-CRM-D (TAB:density): for each density level, the fraction of
%       instances with overlapping regions, the feasibility success rate of the
%       NAIVE single-region offset versus the MERGE-aware offset (which treats an
%       overlapping cluster as one component of D = union_l D_l), the number of
%       merge offset updates k, and the optimality gap Delta of the feasible point
%       relative to the unconstrained optimum f(r*)=sum_i c_i||p_i-r*||^2.
%
%   (2) Distributed H-DRMi-D (TAB:density_dist): feasibility, mean offset
%       applications, consensus rounds, and the energy gap of the distributed
%       feasible point to the centralized (merge) H-CRM-D solution.
%
%   SELF-CONTAINED: a square region R=[0,L]^2 holds D convex denied polytopes
%   (regular hexagons of circumradius rho) placed at random; D is chosen so the
%   NOMINAL coverage equals each 'coverage' level, so overlap grows with density.
%   Feasibility, intersection and clustering use MATLAB polyshape.
%
%   Defaults: coverage = [5 15 30 50] (%), nC = 300, nD = 150.
%   opts: N (30), L (20), rho (2.5), commRange (7), seed (4242), kmax (60),
%         step (0.12*rho), tol (1e-3), maxIter (3000), reRounds (5).
%
%   NOTE: this reproduces the prototype's TRENDS; absolute values depend on the
%   knobs above (rho, step, commRange, kmax). If a column lands off the target,
%   calibrate: lower commRange or raise kmax to drop naive success faster; scale
%   step/rho to move k and the gaps.
%
%   See also EXPERIMENTGAPSWEEP, EXPERIMENTBENCHMARK.

    arguments
        coverage (1,:) double = [5 15 30 50]
        nC       (1,1) double = 300
        nD       (1,1) double = 300
        opts.N         (1,1) double = 30
        opts.L         (1,1) double = 20
        opts.rho       (1,1) double = 2.5
        opts.commRange (1,1) double = 7
        opts.seed      (1,1) double = 4242
        opts.kmax      (1,1) double = 60
        opts.step      (1,1) double = 0
        opts.tol       (1,1) double = 1e-3
        opts.maxIter   (1,1) double = 3000
        opts.reRounds  (1,1) double = 5
    end
    if opts.step == 0, opts.step = 0.12 * opts.rho; end

    L = opts.L;  rho = opts.rho;  N = opts.N;
    areaHex = (3*sqrt(3)/2) * rho^2;
    nL = numel(coverage);
    z  = nan(nL,1);

    % centralized table
    C.cov = coverage(:);  C.overlap = z;  C.succNaive = z;  C.succMerge = z;  C.k = z;  C.dE = z;
    % distributed table
    D.cov = coverage(:);  D.overlap = z;  D.feasible = z;  D.offsets = z;  D.rounds = z;  D.gap = z;

    for li = 1:nL
        Dcount = max(1, round(coverage(li)/100 * L^2 / areaHex));

        % ---------------- centralized: nC instances ----------------
        ov = nan(nC,1); sn = nan(nC,1); sm = nan(nC,1); kk = nan(nC,1); de = nan(nC,1);
        for s = 1:nC
            rng(opts.seed + 100000*li + s);
            [hex, U, anyOv] = makeField(L, rho, Dcount);
            [P, w] = makeAgents(N, L);
            rStar = sum(w.*P,1) ./ sum(w);

            ov(s) = anyOv;
            [~, ~, okN]      = naiveOffset(rStar, hex, U, L, opts);
            [rM, kM, okM]    = mergeOffset(rStar, hex, U, L, opts);
            sn(s) = okN;  sm(s) = okM;  kk(s) = kM;
            de(s) = 100 * (fObj(rM,P,w) - fObj(rStar,P,w)) / max(fObj(rStar,P,w), eps);
        end
        C.overlap(li)   = 100*mean(ov,'omitnan');
        C.succNaive(li) = 100*mean(sn,'omitnan');
        C.succMerge(li) = 100*mean(sm,'omitnan');
        C.k(li)         = mean(kk,'omitnan');
        C.dE(li)        = mean(de(sm==1),'omitnan');     % gap over feasible merge points

        % ---------------- distributed: nD instances ----------------
        ovd = nan(nD,1); fe = nan(nD,1); of = nan(nD,1); rd = nan(nD,1); gp = nan(nD,1);
        for s = 1:nD
            rng(opts.seed + 200000*li + s);
            [hex, U, anyOv] = makeField(L, rho, Dcount);
            [P, w] = makeAgents(N, L);
            A = commGraph(P, opts.commRange);

            % centralized merge reference on the SAME instance
            rStar = sum(w.*P,1) ./ sum(w);
            rRef  = mergeOffset(rStar, hex, U, L, opts);

            % local proposals = neighborhood weighted centroids
            R0 = localProposals(P, w, A);
            [rbar, it0] = consensusAgree(R0, w, A, opts.tol, opts.maxIter);

            % check - offset - reconsensus
            r = rbar;  nOff = 0;  rounds = it0;
            while isinterior(U, r(1), r(2)) && nOff < opts.kmax
                l  = firstContaining(r, hex);
                if isempty(l), break; end
                cl = clusterCentroid(l, hex);
                dir = r - cl;  if norm(dir) < 1e-9, dir = randn(1,2); end
                r = r + opts.step * dir/norm(dir);
                nOff = nOff + 1;  rounds = rounds + opts.reRounds;   % agents re-agree
            end
            ovd(s) = anyOv;
            fe(s)  = double(~isinterior(U, r(1), r(2)) && inBox(r, L));
            of(s)  = nOff;  rd(s) = rounds;
            gp(s)  = 100 * (fObj(r,P,w) - fObj(rRef,P,w)) / max(fObj(rRef,P,w), eps);
        end
        D.overlap(li)  = 100*mean(ovd,'omitnan');
        D.feasible(li) = 100*mean(fe,'omitnan');
        D.offsets(li)  = mean(of,'omitnan');
        D.rounds(li)   = mean(rd,'omitnan');
        D.gap(li)      = mean(gp,'omitnan');

        fprintf(['cov=%2d%% (D=%d): overlap=%.0f%% | naive=%.1f%% merge=%.1f%% k=%.1f dE=%.2f%% ' ...
                 '|| feas=%.0f%% off=%.1f rounds=%.0f gap=%.2f%%\n'], ...
            coverage(li), Dcount, C.overlap(li), C.succNaive(li), C.succMerge(li), C.k(li), C.dE(li), ...
            D.feasible(li), D.offsets(li), D.rounds(li), D.gap(li));
    end

    out.C = C;  out.D = D;  out.coverage = coverage(:);
    printTables(C, D, nC, nD);
    printLatex(C, D);
end

% =====================================================================
function [hex, U, anyOv] = makeField(L, rho, D)
%MAKEFIELD  D random regular hexagons (circumradius rho) in [0,L]^2.
    ang = (0:5)*(pi/3) + pi/6;                         % vertex angles
    vx  = rho*cos(ang);  vy = rho*sin(ang);
    hex = struct('center',{},'poly',{});
    for d = 1:D
        c = rho + (L-2*rho).*rand(1,2);
        hex(d).center = c; %#ok<AGROW>
        hex(d).poly   = polyshape(c(1)+vx, c(2)+vy); %#ok<AGROW>
    end
    U = hex(1).poly;
    for d = 2:D, U = union(U, hex(d).poly); end
    % any overlapping pair?
    anyOv = false;
    for a = 1:D-1
        for b = a+1:D
            if area(intersect(hex(a).poly, hex(b).poly)) > 1e-9
                anyOv = true;  break;
            end
        end
        if anyOv, break; end
    end
end

% ---------------------------------------------------------------------
function [P, w] = makeAgents(N, L)
%MAKEAGENTS  Random positions and positive rendezvous weights c_i.
    P = L * rand(N, 2);
    w = 0.5 + 1.5*rand(N, 1);                          % c_i in [0.5, 2]
end

% ---------------------------------------------------------------------
function v = fObj(r, P, w)
%FOBJ  Energy-weighted dispersion f(r) = sum_i c_i ||p_i - r||^2.
    v = sum(w .* sum((P - r).^2, 2));
end

% ---------------------------------------------------------------------
function l = firstContaining(r, hex)
%FIRSTCONTAINING  Index of the first hexagon that contains point r, else [].
    l = [];
    for d = 1:numel(hex)
        if isinterior(hex(d).poly, r(1), r(2)), l = d;  return; end
    end
end

% ---------------------------------------------------------------------
function comp = componentOf(l, hex)
%COMPONENTOF  Overlap-connected cluster (indices) containing region l.
    D = numel(hex);
    A = false(D);
    for a = 1:D-1
        for b = a+1:D
            if area(intersect(hex(a).poly, hex(b).poly)) > 1e-9
                A(a,b) = true;  A(b,a) = true;
            end
        end
    end
    seen = false(1,D);  stack = l;  seen(l) = true;
    while ~isempty(stack)
        x = stack(end);  stack(end) = [];
        nb = find(A(x,:) & ~seen);
        seen(nb) = true;  stack = [stack, nb]; %#ok<AGROW>
    end
    comp = find(seen);
end

% ---------------------------------------------------------------------
function cl = clusterCentroid(l, hex)
%CLUSTERCENTROID  Centroid of the overlap cluster containing region l.
    comp = componentOf(l, hex);
    cl = mean(vertcat(hex(comp).center), 1);
end

% ---------------------------------------------------------------------
function [r, k, ok] = naiveOffset(r0, hex, U, L, opts)
%NAIVEOFFSET  Single-region offset (one-shot): push out of the ONE containing
%   region only, then stop. In an overlapping cluster the escaped point can land
%   in a neighbor, so feasibility (checked against the whole union) then fails --
%   this is the failure mode the merge-aware variant removes.
    r = r0;  k = 0;
    if ~isinterior(U, r(1), r(2)), ok = inBox(r, L);  return; end     % no offset needed
    l = firstContaining(r, hex);
    dir = r - hex(l).center;  if norm(dir) < 1e-9, dir = randn(1,2); end
    dir = dir / norm(dir);
    subMax = ceil(3*opts.rho/opts.step);  sub = 0;
    while isinterior(hex(l).poly, r(1), r(2)) && sub < subMax
        r = r + opts.step * dir;  sub = sub + 1;
    end
    k = 1;
    ok = ~isinterior(U, r(1), r(2)) && inBox(r, L);                    % fails if now in a neighbor
end

% ---------------------------------------------------------------------
function [r, k, ok] = mergeOffset(r0, hex, U, L, opts)
%MERGEOFFSET  Merge-aware offset: step outward from the whole overlap cluster.
    r = r0;  k = 0;
    if ~isinterior(U, r(1), r(2)), ok = inBox(r, L);  return; end
    l  = firstContaining(r, hex);
    cl = clusterCentroid(l, hex);
    dir = r0 - cl;  if norm(dir) < 1e-9, dir = randn(1,2); end
    dir = dir / norm(dir);
    while isinterior(U, r(1), r(2)) && k < opts.kmax
        r = r + opts.step * dir;  k = k + 1;
    end
    if isinterior(U, r(1), r(2)) || ~inBox(r, L)       % radial blocked: try opposite
        r = r0;  dir = -dir;  k2 = 0;
        while isinterior(U, r(1), r(2)) && k2 < opts.kmax
            r = r + opts.step * dir;  k2 = k2 + 1;
        end
        k = k + k2;
    end
    ok = ~isinterior(U, r(1), r(2)) && inBox(r, L);
end

% ---------------------------------------------------------------------
function tf = inBox(r, L)
    tf = all(r >= 0) && all(r <= L);
end

% ---------------------------------------------------------------------
function A = commGraph(P, range)
%COMMGRAPH  Undirected range-limited communication graph (no self-loops).
    N = size(P,1);
    sq = sum(P.^2, 2);
    Dd = sqrt(max(sq + sq.' - 2*(P*P.'), 0));
    Dd(1:N+1:end) = Inf;
    A = Dd <= range;
end

% ---------------------------------------------------------------------
function R0 = localProposals(P, w, A)
%LOCALPROPOSALS  Each agent's neighborhood weighted centroid (incl. self).
    N = size(P,1);  R0 = P;
    for i = 1:N
        nb = [find(A(i,:)), i];
        R0(i,:) = sum(w(nb).*P(nb,:), 1) ./ sum(w(nb));
    end
end

% ---------------------------------------------------------------------
function [rbar, iters] = consensusAgree(R0, w, A, tol, maxIter)
%CONSENSUSAGREE  Self-normalizing power-weighted consensus to a common point.
    R = R0;  N = size(R0,1);  iters = maxIter;
    for it = 1:maxIter
        Rn = R;
        for i = 1:N
            nb = find(A(i,:));
            if isempty(nb), continue; end
            wn = w(nb).';  denom = w(i) + sum(wn);
            Rn(i,:) = R(i,:) + (wn*(R(nb,:) - R(i,:))) / denom;
        end
        d = max(vecnorm(Rn - R, 2, 2));
        R = Rn;
        if d < tol, iters = it;  break; end
    end
    rbar = mean(R, 1);
end

% ---------------------------------------------------------------------
function printTables(C, D, nC, nD)
    fprintf('\n=== H-CRM-D (centralized), %d instances/level ===\n', nC);
    fprintf('%-8s %-9s %-11s %-11s %-6s %-8s\n', ...
        'cov(%)','ovlp(%)','naive(%)','merge(%)','k','dE(%)');
    for i = 1:numel(C.cov)
        fprintf('%-8d %-9.0f %-11.1f %-11.1f %-6.1f %-8.2f\n', ...
            C.cov(i), C.overlap(i), C.succNaive(i), C.succMerge(i), C.k(i), C.dE(i));
    end
    fprintf('\n=== H-DRMi-D (distributed), %d instances/level ===\n', nD);
    fprintf('%-8s %-9s %-11s %-9s %-8s %-9s\n', ...
        'cov(%)','ovlp(%)','feas(%)','offsets','rounds','gap(%)');
    for i = 1:numel(D.cov)
        fprintf('%-8d %-9.0f %-11.0f %-9.1f %-8.0f %-9.2f\n', ...
            D.cov(i), D.overlap(i), D.feasible(i), D.offsets(i), D.rounds(i), D.gap(i));
    end
    fprintf('================================================\n');
end

% ---------------------------------------------------------------------
function printLatex(C, D)
    fprintf('\n%% --- LaTeX rows for tab:density (H-CRM-D) ---\n');
    for i = 1:numel(C.cov)
        fprintf('$%d$ & $%.0f$ & $%.1f$ & $%.0f$ & $%.1f$ & $%.2f$\\\\\n', ...
            C.cov(i), C.overlap(i), C.succNaive(i), C.succMerge(i), C.k(i), C.dE(i));
    end
    fprintf('\n%% --- LaTeX rows for tab:density_dist (H-DRMi-D) ---\n');
    for i = 1:numel(D.cov)
        fprintf('$%d$ & $%.0f$ & $%.0f$ & $%.1f$ & $%.0f$ & $%.2f$\\\\\n', ...
            D.cov(i), D.overlap(i), D.feasible(i), D.offsets(i), D.rounds(i), D.gap(i));
    end
    fprintf('%% ------------------------------------------------\n');
end
