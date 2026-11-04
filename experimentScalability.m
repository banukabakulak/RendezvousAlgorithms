function out = experimentScalability(Nvalues, nInst, baseParams, N0)
%EXPERIMENTSCALABILITY  C-MAS vs D-MAS scalability over team size N.
%
%   out = EXPERIMENTSCALABILITY(Nvalues, nInst, baseParams, N0) sweeps the number
%   of agents N and, for each N, reports the four scalability metrics plotted in
%   Fig.~\ref{fig:ScalabilityPlot}:
%
%     (a) per-cycle algorithm runtime (ms): the wall-clock cost of ONE rendezvous
%         decision. For C-MAS this is the full H-CRM-D solve over all N agents; for
%         D-MAS it is the per-agent critical path max_i (local H-CRM-D solve over
%         agent i's neighborhood), i.e. the parallel per-cycle bottleneck
%         O(max_i (M^(i))^4 + n D^(i)) rather than the centralized O(M^4 + nD).
%     (b) convergence time Tc (s);
%     (c) communication overhead (Kbits);
%     (d) mission energy (mAh).
%
%   Metrics (b)-(d) are taken from a full simulation of each controller via
%   RUNCMAS/RUNDMAS and COMPUTEMETRICS/COMPUTEMETRICSDMAS. To keep the comparison
%   fair the instance is regenerated (GENERATEINSTANCE resets rng(seed)) before
%   each controller so both run from an identical RNG state.
%
%   The cubic arena is scaled with N (side proportional to N^{1/3}, with the RoI
%   and denied region scaled to match) so the mean communication degree stays
%   roughly constant and the multi-hop regime of Sec.~\ref{sec:InstanceGeneration}
%   is preserved as the team grows; the communication ranges are held fixed.
%
%   Defaults: Nvalues = [10 25 50 100 200], nInst = 10, baseParams =
%   defaultParams(), N0 = 50 (the reference size at which the base arena is used).
%
%   OUT has column vectors over Nvalues: N, degMean, tCMAS, tDMAS (ms), TcC, TcD
%   (s), ovC, ovD (Kbits), eC, eD (mAh), and OUT.raw with per-instance values.
%
%   NOTE on metric field names: COMPUTEMETRICS / COMPUTEMETRICSDMAS field names
%   are resolved defensively (see PICK below). If a metric comes back NaN, add the
%   actual field name to the candidate list in PICK.
%
%   See also RUNCMAS, RUNDMAS, COMPUTEMETRICS, COMPUTEMETRICSDMAS, PLOTSCALABILITY.

    arguments
        Nvalues    (1,:) double = [10 25 50 100 200]
        nInst      (1,1) double = 10
        baseParams (1,1) struct = defaultParams()
        N0         (1,1) double = 50
    end

    baseParams.controller.logGraph = false;            % not needed; saves memory
    W0 = baseParams.arena.width;  H0 = baseParams.arena.height;
    roi0 = baseParams.region.roi; denied0 = baseParams.region.denied;

    nN = numel(Nvalues);
    z  = nan(nN,1);
    out.N = Nvalues(:);  out.degMean = z;
    out.tCMAS = z; out.tDMAS = z;  out.TcC = z; out.TcD = z;
    out.ovC   = z; out.ovD   = z;  out.eC  = z; out.eD  = z;
    out.raw = struct('N',{},'deg',{},'tCMAS',{},'tDMAS',{}, ...
                     'TcC',{},'TcD',{},'ovC',{},'ovD',{},'eC',{},'eD',{});

    for ni = 1:nN
        N  = Nvalues(ni);
        sc = (N / N0)^(1/3);                           % cubic-arena side scaling

        dg = nan(nInst,1);
        tC = nan(nInst,1); tD = nan(nInst,1);
        TcC = nan(nInst,1); TcD = nan(nInst,1);
        ovC = nan(nInst,1); ovD = nan(nInst,1);
        eC  = nan(nInst,1); eD  = nan(nInst,1);

        for s = 1:nInst
            p = baseParams;
            p.sim.N         = N;
            p.arena.width   = W0 * sc;
            p.arena.height  = H0 * sc;
            p.region.roi    = roi0    * sc;
            p.region.denied = denied0 * sc;
            p.sim.seed      = 70000 + 1000*ni + s;     % reproducible

            % ----- instance for the per-cycle runtime micro-benchmark -----
            inst = generateInstance(p);
            [A, deg] = localAdjacency(inst);
            dg(s) = mean(deg);

            % (a) C-MAS per-cycle: full H-CRM-D solve over all N agents
            tC(s) = timeCall(@() solveRendezvousCentral(p, inst), 0.01) * 1e3;

            % (a) D-MAS per-cycle: critical path = max over agents of the local
            %     H-CRM-D solve on that agent's neighborhood (incl. self)
            tLoc = nan(N,1);
            for i = 1:N
                nb  = [find(A(i,:)), i];
                sub = subInstance(inst, nb);
                tLoc(i) = timeCall(@() solveRendezvousCentral(p, sub), 0.004) * 1e3;
            end
            tD(s) = max(tLoc);

            % ----- full-simulation metrics (fair RNG: regenerate each run) -----
            instC = generateInstance(p);
            resC  = runCMAS(p, instC);
            mC    = computeMetrics(p, resC);

            instD = generateInstance(p);
            resD  = runDMAS(p, instD);
            mD    = computeMetricsDMAS(p, resD);

            TcC(s) = pick(mC, {'avgConsensusTime','Tc','convergenceTime'});
            TcD(s) = pick(mD, {'avgConsensusTime','Tc','convergenceTime'});
            ovC(s) = pick(mC, {'totalOverhead','overhead','O'});
            ovD(s) = pick(mD, {'totalOverhead','overhead','O'});
            eC(s)  = pick(mC, {'totalEnergy','energy','E'});
            eD(s)  = pick(mD, {'totalEnergy','energy','E'});

            out.raw(end+1) = struct('N',N,'deg',dg(s),'tCMAS',tC(s),'tDMAS',tD(s), ...
                'TcC',TcC(s),'TcD',TcD(s),'ovC',ovC(s),'ovD',ovD(s), ...
                'eC',eC(s),'eD',eD(s)); %#ok<AGROW>
        end

        out.degMean(ni) = mean(dg,'omitnan');
        out.tCMAS(ni)   = mean(tC,'omitnan');   out.tDMAS(ni) = mean(tD,'omitnan');
        out.TcC(ni)     = mean(TcC,'omitnan');  out.TcD(ni)   = mean(TcD,'omitnan');
        out.ovC(ni)     = mean(ovC,'omitnan');  out.ovD(ni)   = mean(ovD,'omitnan');
        out.eC(ni)      = mean(eC,'omitnan');   out.eD(ni)    = mean(eD,'omitnan');

        fprintf(['N=%3d (deg=%.1f): runtime C=%.3f D=%.3f ms | Tc C=%.1f D=%.1f s | ' ...
                 'O C=%.0f D=%.0f Kbits | E C=%.0f D=%.0f mAh\n'], ...
            N, out.degMean(ni), out.tCMAS(ni), out.tDMAS(ni), out.TcC(ni), out.TcD(ni), ...
            out.ovC(ni), out.ovD(ni), out.eC(ni), out.eD(ni));
    end

    printTable(out);
end

% =====================================================================
function [A, deg] = localAdjacency(inst)
%LOCALADJACENCY  Mutual (undirected) comm graph from positions and ranges.
%   Link i-j iff dist(i,j) <= min(range_i, range_j). Matches GENERATEINSTANCE.
    P = inst.pos;  r = inst.cap.commRange;  N = size(P,1);
    sq = sum(P.^2,2);
    D  = sqrt(max(sq + sq.' - 2*(P*P.'), 0));
    D(1:N+1:end) = Inf;
    A  = D <= min(r, r.');
    deg = sum(A, 2);
end

% ---------------------------------------------------------------------
function sub = subInstance(inst, idx)
%SUBINSTANCE  Restrict an instance to the agent index set IDX (keeps regions).
    idx = unique(idx(:));
    sub = inst;
    sub.agentType = inst.agentType(idx);
    sub.id        = inst.id(idx);
    sub.pos       = inst.pos(idx,:);
    sub.vel       = inst.vel(idx,:);
    sub.speed     = inst.speed(idx);
    sub.power     = inst.power(idx);
    sub.weight    = inst.weight(idx);
    for f = ["powerA","powerB","vMin","vMax","aMin","aMax","commRange","battery"]
        sub.cap.(f) = inst.cap.(f)(idx);
    end
    for f = ["voltage","qDischarged","soc"]
        if isfield(inst.battery, f), sub.battery.(f) = inst.battery.(f)(idx); end
    end
    w = sub.weight;
    sub.centroid = sum(w .* sub.pos, 1) ./ sum(w);
    % sub.roi and sub.denied are inherited unchanged.
end

% ---------------------------------------------------------------------
function t = timeCall(fn, minBudget)
%TIMECALL  Stable per-call time (s): accumulate calls until MINBUDGET seconds.
    fn();                                              % warm up
    k = 1;  el = 0;
    while el < minBudget && k <= 1e6
        t0 = tic;
        for i = 1:k, fn(); end
        el = toc(t0);
        if el < minBudget, k = k * 4; end
    end
    t = el / k;
end

% ---------------------------------------------------------------------
function v = pick(s, names)
%PICK  First existing field of struct S among NAMES, else NaN (with one warning).
    for k = 1:numel(names)
        if isfield(s, names{k}), v = s.(names{k});  return; end
    end
    v = NaN;
    warning('experimentScalability:field', ...
            'None of {%s} found in metrics struct; reporting NaN.', strjoin(names, ', '));
end

% ---------------------------------------------------------------------
function printTable(out)
    fprintf('\n=== Scalability: C-MAS vs D-MAS (mean over instances) ===\n');
    fprintf('%-5s %-6s %-9s %-9s %-8s %-8s %-10s %-10s %-9s %-9s\n', ...
        'N','deg','tC(ms)','tD(ms)','TcC(s)','TcD(s)','OC(Kb)','OD(Kb)','EC(mAh)','ED(mAh)');
    for i = 1:numel(out.N)
        fprintf('%-5d %-6.1f %-9.3f %-9.3f %-8.1f %-8.1f %-10.0f %-10.0f %-9.0f %-9.0f\n', ...
            out.N(i), out.degMean(i), out.tCMAS(i), out.tDMAS(i), out.TcC(i), out.TcD(i), ...
            out.ovC(i), out.ovD(i), out.eC(i), out.eD(i));
    end
    fprintf('=========================================================\n');
    fprintf('overhead ratio D/C at N=%d: %.1fx\n', out.N(end), out.ovD(end)/max(out.ovC(end),eps));
end
