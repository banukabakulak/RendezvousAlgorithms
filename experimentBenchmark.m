function out = experimentBenchmark(Mvalues, nInst, opts)
%EXPERIMENTBENCHMARK  Runtime benchmark of CRM-OPT versus general QP solvers.
%
%   out = EXPERIMENTBENCHMARK(Mvalues, nInst) benchmarks the structure-exploiting
%   CRM-OPT rendezvous solver against MATLAB's general quadratic-programming
%   solvers on constrained CRM instances of increasing size M (the number of
%   region facets / linear constraints), over nInst random instances per size.
%
%   Because the CRM objective has the isotropic Hessian 2(sum_i c_i) I, the
%   constrained CRM is the Euclidean projection of the weighted centroid r* onto
%   the region R. Each instance is therefore a projection problem:
%
%       min_r  (sum_i c_i) || r - r* ||^2      s.t.   A r <= b,
%
%   with R a convex polygon of M facets (a near-regular M-gon of unit circumradius
%   centred at the origin) and r* drawn OUTSIDE R so the projection is active (the
%   interesting, constraint-binding regime). The three solvers compared are:
%
%     CRM-OPT       O(M) projection of r* onto the convex polygon (this file);
%                   the structure-exploiting kernel used by SOLVERENDEZVOUSCENTRAL.
%     quadprog-AS   active-set QP (quadprog 'active-set'; if that algorithm is not
%                   available in this MATLAB release, lsqlin 'active-set', which
%                   solves the identical constrained least-squares projection).
%     quadprog-IP   interior-point-convex QP (quadprog 'interior-point-convex'),
%                   the general-purpose solver that ignores the isotropic structure.
%
%   An optional commercial solver (Gurobi / MOSEK) is timed when a solver handle is
%   supplied via the 'Commercial' option (see below); otherwise that column is
%   omitted.
%
%   All solvers are checked to return the same optimum: OUT.gap is the maximum
%   point deviation (m) of any QP solver from the CRM-OPT point over all instances,
%   which confirms exactness.
%
%   Name-Value options:
%     'Reps'        min accumulated time (s) per timing (default 0.02) -- each
%                   solver is called repeatedly until this budget is reached and
%                   the per-call time is reported, so sub-millisecond kernels are
%                   measured reliably.
%     'Commercial'  function handle r = f(H,f,A,b,r0) for a commercial solver;
%                   default [] (column omitted).
%     'ComName'     column label for the commercial solver (default 'Gurobi').
%     'Seed'        base RNG seed (default 12345), for reproducibility.
%
%   Defaults: Mvalues = [4 8 16 32 64 96], nInst = 50.
%
%   OUT fields (column vectors over Mvalues): M, tCRM, tAS, tIP, tCom (ms, mean),
%   gap (m, mean max-deviation), nInst; OUT.raw holds per-instance values and
%   OUT.comName / OUT.hasCom describe the commercial column. A readable table and
%   paste-ready LaTeX rows are printed at the end.
%
%   See also SOLVERENDEZVOUSCENTRAL, PLOTBENCHMARK.

    arguments
        Mvalues (1,:) double = [4 8 16 32 64 96]
        nInst   (1,1) double = 50
        opts.Reps       (1,1) double        = 0.02
        opts.Commercial                     = []
        opts.ComName    (1,:) char          = 'Gurobi'
        opts.Seed       (1,1) double        = 12345
    end

    hasCom = ~isempty(opts.Commercial);
    nM = numel(Mvalues);
    z  = nan(nM,1);
    out.M = Mvalues(:);
    out.tCRM = z; out.tAS = z; out.tIP = z; out.tCom = z; out.gap = z;
    out.nInst = nInst * ones(nM,1);
    out.hasCom = hasCom; out.comName = opts.ComName;
    out.raw = struct('M',{},'tCRM',{},'tAS',{},'tIP',{},'tCom',{},'gap',{});

    asBackend = '';                                   % resolved on first use

    for mi = 1:nM
        M = Mvalues(mi);
        tC = nan(nInst,1); tA = nan(nInst,1); tI = nan(nInst,1);
        tK = nan(nInst,1); gp = nan(nInst,1);

        for s = 1:nInst
            rng(opts.Seed + 1000*mi + s);             % reproducible per instance

            % ---- region: near-regular M-gon, unit circumradius, centred at 0 ----
            th = 2*pi*((0:M-1).')/M + 0.15*(2*pi/M)*(rand(M,1)-0.5);  % slight jitter
            th = sort(th);
            V  = [cos(th), sin(th)];                  % M x 2 vertices (CCW)

            % halfspaces A r <= b from the edges (outward normals)
            [A, b] = polygonHalfspaces(V);

            % ---- weighted centroid r* placed OUTSIDE the region ----
            u     = randn(1,2); u = u / norm(u);
            rStar = (1.5 + 1.5*rand) * u;             % radius in [1.5, 3] > 1
            rStar = rStar(:);

            S  = 0.5 + rand;                          % sum_i c_i > 0 (scales objective only)
            H  = 2*S*eye(2);
            f  = -2*S*rStar;
            r0 = zeros(2,1);                          % origin is feasible (interior)

            % ---- CRM-OPT: O(M) projection onto the convex polygon ----
            fnC    = @() crmProject(rStar, V);
            rC     = fnC();
            tC(s)  = timeCall(fnC, opts.Reps) * 1e3;  % ms

            % ---- active-set QP (quadprog 'active-set', else lsqlin 'active-set') ----
            [fnA, asBackend] = makeActiveSet(H, f, A, b, r0, rStar, asBackend);
            rA     = fnA();
            tA(s)  = timeCall(fnA, opts.Reps) * 1e3;

            % ---- interior-point-convex QP (general solver) ----
            optIP  = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
            fnI    = @() quadprog(H, f, A, b, [], [], [], [], [], optIP);
            rI     = fnI();
            tI(s)  = timeCall(fnI, opts.Reps) * 1e3;

            % ---- optional commercial solver ----
            dev = [norm(rA(:)-rC(:)), norm(rI(:)-rC(:))];
            if hasCom
                fnK    = @() opts.Commercial(H, f, A, b, r0);
                rK     = fnK();
                tK(s)  = timeCall(fnK, opts.Reps) * 1e3;
                dev(end+1) = norm(rK(:)-rC(:)); %#ok<AGROW>
            end
            gp(s) = max(dev);

            out.raw(end+1) = struct('M',M,'tCRM',tC(s),'tAS',tA(s),'tIP',tI(s), ...
                                    'tCom',tK(s),'gap',gp(s)); %#ok<AGROW>
        end

        out.tCRM(mi) = mean(tC,'omitnan');
        out.tAS(mi)  = mean(tA,'omitnan');
        out.tIP(mi)  = mean(tI,'omitnan');
        out.tCom(mi) = mean(tK,'omitnan');
        out.gap(mi)  = mean(gp,'omitnan');

        fprintf('M=%3d (%d inst): CRM-OPT=%.4f ms  AS=%.3f ms  IP=%.3f ms', ...
            M, nInst, out.tCRM(mi), out.tAS(mi), out.tIP(mi));
        if hasCom, fprintf('  %s=%.3f ms', opts.ComName, out.tCom(mi)); end
        fprintf('  | max gap=%.2e m\n', out.gap(mi));
    end

    out.asBackend = asBackend;
    if ~isempty(asBackend)
        fprintf('\n(active-set column computed with: %s)\n', asBackend);
    end
    printReadable(out);
    printLatex(out);
end

% =====================================================================
function r = crmProject(p, V)
%CRMPROJECT  Euclidean projection of point P onto the convex polygon V (CCW).
%   If P is inside, return P; else return the nearest point on the boundary over
%   all M edges. Fully vectorized (no per-edge loop), so the kernel cost is
%   essentially flat in M -- the structure-exploiting CRM-OPT kernel used by
%   SOLVERENDEZVOUSCENTRAL.
    p = p(:).';
    M = size(V,1);
    W = V([2:M, 1], :);                               % next vertex of each edge
    E = W - V;                                        % edge vectors
    % inside test (CCW): P is left of / on every edge
    cross = E(:,1).*(p(2)-V(:,2)) - E(:,2).*(p(1)-V(:,1));
    if all(cross >= -1e-12)
        r = p(:);  return;
    end
    L2 = sum(E.^2, 2);
    t  = max(0, min(1, sum((p - V).*E, 2) ./ max(L2, eps)));
    Q  = V + t.*E;                                    % foot of P on each edge
    d  = sum((p - Q).^2, 2);
    [~, k] = min(d);
    r  = Q(k, :).';
end

% ---------------------------------------------------------------------
function [A, b] = polygonHalfspaces(V)
%POLYGONHALFSPACES  Outward-normal halfspace form A r <= b of a CCW polygon.
    M = size(V,1);  A = zeros(M,2);  b = zeros(M,1);
    c = mean(V,1);                                    % interior point
    for k = 1:M
        a = V(k,:);  bnd = V(mod(k,M)+1,:);  e = bnd - a;
        n = [e(2), -e(1)];  n = n / norm(n);          % edge normal
        if n*(c-a).' > 0, n = -n; end                 % force OUTWARD
        A(k,:) = n;  b(k) = n*a.';
    end
end

% ---------------------------------------------------------------------
function [fn, backend] = makeActiveSet(H, f, A, b, r0, rStar, backend)
%MAKEACTIVESET  Return a timing handle for the active-set solver, resolving the
%   backend once: quadprog 'active-set' if this release supports it, else lsqlin
%   'active-set' on the equivalent least-squares projection.
    n = numel(r0);
    if isempty(backend)
        try
            o = optimoptions('quadprog','Algorithm','active-set','Display','off');
            quadprog(H, f, A, b, [], [], [], [], r0, o);     % probe
            backend = 'quadprog active-set';
        catch
            backend = 'lsqlin active-set';
        end
    end
    if strcmp(backend, 'quadprog active-set')
        o  = optimoptions('quadprog','Algorithm','active-set','Display','off');
        fn = @() quadprog(H, f, A, b, [], [], [], [], r0, o);
    else
        o  = optimoptions('lsqlin','Algorithm','active-set','Display','off');
        C  = eye(n);  d = rStar(:);
        fn = @() lsqlin(C, d, A, b, [], [], [], [], r0, o);
    end
end

% ---------------------------------------------------------------------
function t = timeCall(fn, minBudget)
%TIMECALL  Stable per-call time (s): accumulate calls until MINBUDGET seconds.
    fn();                                             % warm up (JIT)
    k = 1;  el = 0;
    while el < minBudget && k <= 1e7
        t0 = tic;
        for i = 1:k, fn(); end
        el = toc(t0);
        if el < minBudget, k = k * 4; end
    end
    t = el / k;
end

% ---------------------------------------------------------------------
function printReadable(out)
    fprintf('\n=== CRM-OPT runtime benchmark (mean over instances) ===\n');
    if out.hasCom
        fprintf('%-5s %-11s %-12s %-12s %-12s %-10s\n', ...
            'M','CRM-OPT(ms)','AS(ms)','IP(ms)',[out.comName '(ms)'],'gap(m)');
        for i = 1:numel(out.M)
            fprintf('%-5d %-11.4f %-12.3f %-12.3f %-12.3f %-10.1e\n', ...
                out.M(i), out.tCRM(i), out.tAS(i), out.tIP(i), out.tCom(i), out.gap(i));
        end
    else
        fprintf('%-5s %-11s %-12s %-12s %-10s\n', ...
            'M','CRM-OPT(ms)','AS(ms)','IP(ms)','gap(m)');
        for i = 1:numel(out.M)
            fprintf('%-5d %-11.4f %-12.3f %-12.3f %-10.1e\n', ...
                out.M(i), out.tCRM(i), out.tAS(i), out.tIP(i), out.gap(i));
        end
    end
    fprintf('=======================================================\n');
end

% ---------------------------------------------------------------------
function printLatex(out)
%PRINTLATEX  Paste-ready LaTeX body rows for tab:benchmark.
    fprintf('\n%% --- LaTeX rows for tab:benchmark (max gap = %.1e m) ---\n', max(out.gap));
    for i = 1:numel(out.M)
        if out.hasCom
            fprintf('$%d$ & $%.4f$ & $%.3f$ & $%.3f$ & $%.3f$ \\\\\n', ...
                out.M(i), out.tCRM(i), out.tAS(i), out.tIP(i), out.tCom(i));
        else
            fprintf('$%d$ & $%.4f$ & $%.3f$ & $%.3f$ \\\\\n', ...
                out.M(i), out.tCRM(i), out.tAS(i), out.tIP(i));
        end
    end
    fprintf('%% ----------------------------------------------------\n');
end
